#!/usr/bin/env python3

from Bio import Entrez
from box import Box
from datetime import datetime
from dateutil.parser import parse as parse_date
from pprint import pprint
from time import sleep
from typing import List
from pytz import timezone
import json
import os
import random
import time
import twitter
import xmltodict
import json


TIMEZONE = timezone('EST')

def get_nodes(start_time: datetime) -> List[Box]:
    """
    Returns a list of nodes that have been created/updated/published after a
    specified date
    """

    start_time = start_time.strftime("%Y/%m/%d")

    # I haven't been able to figure out how to query by hour so we have to
    # parse the times and filter manually. No biggie.
    term = f'"{start_time}"[EDAT] : "3000"[EDAT]'

    node_list = Box(
        xmltodict.parse(
            Entrez.esearch("taxonomy", term, retmax=100, email="austin@onecodex.com").read()
        )
    )

    with open('node-list.json', 'w') as handle:
        json.dump(node_list.to_dict(), handle)

    if node_list.eSearchResult.IdList is None:
        return []
    else:
        tax_ids = node_list.eSearchResult.IdList.Id
        response_data = \
            xmltodict.parse(
                Entrez.efetch(
                    db="taxonomy", id=",".join(tax_ids), email="austin@onecodex.com"
                ).read()
            )

        nodes = Box(response_data)

        with open('nodes.json', 'w') as handle:
            json.dump(nodes.to_dict(), handle)

        return [
            Box(
                {
                    "id": node.TaxId,
                    "name": node.ScientificName,
                    "rank": node.Rank,
                    "created_at": parse_date(node.CreateDate).replace(tzinfo=TIMEZONE),
                    "updated_at": parse_date(node.UpdateDate).replace(tzinfo=TIMEZONE),
                    "published_at": parse_date(node.PubDate).replace(tzinfo=TIMEZONE),
                    "lineage": node.Lineage,
                }
            )
            for node in nodes.TaxaSet.Taxon
        ]


def format_tweet_for_node(node, last_time) -> str:

    if node.created_at >= last_time:
        action = "New"
    elif node.published_at >= last_time:
        action = "New"
    elif node.updated_at >= last_time:
        action = "Updated"
    url = f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={node.id}"

    if "Bacteria" in node.lineage:
        emoji = "🦠"
    elif "Fungi" in node.lineage:
        emoji = "🍄"
    elif "Aves" in node.lineage:
        emoji = "🐦"
    elif "Aranea" in node.lineage:
        emoji = "🕷"
    elif "Hymenoptera" in node.lineage:
        emoji = "🐝"
    elif "Formicidae" in node.lineage:
        emoji = "🐜"
    elif "Hexapoda" in node.lineage:
        emoji = "🐞"
    elif "Lepidoptera" in node.lineage:
        emoji = "🦋"
    elif "Viridiplantae" in node.lineage:
        emoji = "🌱"
    elif "Mammalia" in node.lineage and "Rodentia" in node.lineage:
        emoji = "🐁"
    else:
        emoji = None

    if emoji is not None:
        action = f"{action} {emoji}"

    return "\n".join([f"{action}! {node.name} ({node.rank})", url])


def get_twitter():
    api = twitter.Api(
        consumer_key=os.getenv("TWITTER_CONSUMER_KEY"),
        consumer_secret=os.getenv("TWITTER_CONSUMER_SECRET"),
        access_token_key=os.getenv("TWITTER_ACCESS_TOKEN_KEY"),
        access_token_secret=os.getenv("TWITTER_ACCESS_TOKEN_SECRET"),
    )

    return api


def send_tweet(tweet_text, dry_run=True):
    api = get_twitter()
    if dry_run:
        print(tweet_text)
    else:
        print(f"tweeting!")
        print(tweet_text)
        print("\n")

        try:
            api.PostUpdate(tweet_text)
        # probably due to tweeting the same thing twice ...
        except twitter.error.TwitterError as e:
            print(e)
            return


def tweet_nodes(nodes, last_time, dry_run=True, delay=60 * 15):
    """
    Tweets updates. Rate limited (max 1 tweet every 5 minutes)
    """

    print(f"tweeting about {len(nodes)} nodes!")

    for node in nodes:
        send_tweet(format_tweet_for_node(node, last_time), dry_run=dry_run)

        if not dry_run:
            time.sleep(delay)


def main():

    with open("last-time.json") as handle:
        last_time = parse_date(json.load(handle)["last-time"])

    # make sure we're logged in before starting anything
    api = get_twitter()
    api.VerifyCredentials()

    print(f"looking up taxa that were created since {last_time}")

    nodes = get_nodes(start_time=last_time)
    print(f"{len(nodes)} nodes fetch from NCBI (since={last_time})")

    new_nodes = []

    for node in nodes:
        if 'sp.' not in node.name:

            if (node.created_at >= last_time) or (node.published_at >= last_time):
                print(f"[keep] {node.name} created={node.created_at} updated={node.updated_at} published={node.published_at}")
                new_nodes.append(node)
            else:
                pprint([node, last_time, node.created_at >= last_time, node.published_at >= last_time, 'sp.' not in node.name])
        else:
            print(f"[skip] {node.name} {node.created_at}/({node.updated_at})")

    print(f"got {len(new_nodes)} new nodes created/updated since {last_time}")

    with open("last-time.json", "w") as handle:
        json.dump({"last-time": datetime.now(TIMEZONE)}, handle, default=str)

    # this will rate limit
    # if something goes wrong, it will not duplicate tweets
    tweet_nodes(new_nodes, last_time, dry_run=False)


if __name__ == "__main__":
    main()
