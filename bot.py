#!/usr/bin/env python3

from Bio import Entrez
from box import Box
from datetime import datetime
from dateutil.parser import parse as parse_date
from pprint import pprint
from time import sleep
from typing import List
import json
import os
import random
import time
import twitter
import xmltodict


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
        xmltodict.parse(Entrez.esearch("taxonomy", term, email="austin@onecodex.com").read())
    )

    pprint(node_list.to_dict())

    if node_list.eSearchResult.IdList is None:
        return []
    else:
        tax_ids = node_list.eSearchResult.IdList.Id

        nodes = Box(
            xmltodict.parse(
                Entrez.efetch(
                    db="taxonomy", id=",".join(tax_ids), email="austin@onecodex.com"
                ).read()
            )
        )

        return [
            Box(
                {
                    "id": node.TaxId,
                    "name": node.ScientificName,
                    "rank": node.Rank,
                    "created_at": parse_date(node.CreateDate),
                    "updated_at": parse_date(node.UpdateDate),
                    "published_at": parse_date(node.PubDate),
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
        print(api.PostUpdate(tweet_text))


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

    if os.path.exists("last-time.json"):
        with open("last-time.json") as handle:
            last_time = parse_date(json.load(handle)["last-time"])
    else:
        # how to get in UTC?
        last_time = datetime.utcnow()

    # make sure we're logged in before starting anything
    api = get_twitter()
    api.VerifyCredentials()

    print(f"looking up taxa that were created since {last_time}")

    nodes = get_nodes(start_time=last_time)
    print(f"{len(nodes)} nodes fetch from NCBI (since={last_time})")

    new_nodes = [n for n in nodes if (n.created_at >= last_time) or (n.published_at >= last_time)]

    print(f"got {len(new_nodes)} new nodes created/updated since {last_time}")

    with open("last-time.json", "w") as handle:
        json.dump({"last-time": datetime.utcnow()}, handle, default=str)

    # this will rate limit
    # if something goes wrong, it will not duplicate tweets
    tweet_nodes(new_nodes, last_time, dry_run=False)


if __name__ == "__main__":
    main()
