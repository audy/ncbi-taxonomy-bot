#!/usr/bin/env python3

from Bio import Entrez
from box import Box
from datetime import datetime
from dateutil.parser import parse as parse_date
from pprint import pprint
from typing import List
from pytz import timezone
import argparse
import json
import os
import random
import time
import twitter
import xmltodict
import json


TIMEZONE = timezone("EST")

def parse_start_time(start_time_string):
    start_time = datetime.fromisoformat(start_time_string)
    assert start_time.tzinfo is not None
    start_time.replace(tzinfo=TIMEZONE)
    print(start_time)
    return start_time

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--start-time', default=None, type=parse_start_time)
    return parser.parse_args()

def get_nodes(start_time: datetime) -> List[Box]:
    """
    Returns a list of nodes that have been created/updated/published after a
    specified date
    """

    # I haven't been able to figure out how to query by hour so we have to
    # parse the times and filter manually. No biggie.
    term = f'"{start_time.strftime("%Y/%m/%d")}"[EDAT] : "3000"[EDAT]'

    node_list = Box(
        xmltodict.parse(
            Entrez.esearch("taxonomy", term, retmax=100, email="austin@onecodex.com").read()
        )
    )

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
        retval = []

        for node in nodes.TaxaSet.Taxon:
            if 'PubDate' in node:
                pub_date = parse_date(node.PubDate).replace(tzinfo=TIMEZONE)
            else:
                # sometimes nodes don't have a PubDate (even though they're public)
                # so just replace it with now
                pub_date = datetime.now(TIMEZONE)
            retval.append(
                Box({
                        "id": node.TaxId,
                        "name": node.ScientificName,
                        "rank": node.Rank,
                        "created_at": parse_date(node.CreateDate).replace(tzinfo=TIMEZONE),
                        "updated_at": parse_date(node.UpdateDate).replace(tzinfo=TIMEZONE),
                        "published_at": pub_date,
                        "lineage": node.Lineage,
                })
            )
        return retval


def format_tweet_for_node(node, start_time) -> str:
    if node.created_at >= start_time:
        action = "New"
    elif node.published_at >= start_time:
        action = "New"
    elif node.updated_at >= start_time:
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


def tweet_nodes(nodes, start_time, dry_run=True, delay=60 * 15):
    """
    Tweets updates. Rate limited (max 1 tweet every 5 minutes)
    """

    print(f"tweeting about {len(nodes)} nodes!")

    for node in nodes:
        send_tweet(format_tweet_for_node(node, start_time), dry_run=dry_run)

        if not dry_run:
            time.sleep(delay)


def main():
    args = parse_args()
    # make sure we're logged in before starting anything
    api = get_twitter()
    api.VerifyCredentials()

    if args.start_time is None:
        start_time = datetime.now(TIMEZONE)
    else:
        start_time = args.start_time

    print(f"start time: {start_time}")

    while True:
        nodes = [ n for n in get_nodes(start_time) if 'sp.' not in n.name ]

        for node in nodes:
            print(str(node.created_at))

        if len(nodes) == 0:
            print(f"no nodes... I sleep (start_time={start_time})")
        else:
            print(f"{len(nodes)} nodes fetch from NCBI (since={start_time})")

            times = []
            for node in nodes:
                times.extend([node.created_at, node.updated_at, node.published_at])

            tweetable_nodes = [
                n for n in nodes
                if (n.created_at >= start_time) or (n.published_at >= start_time) or (n.updated_at >= start_time)
            ]

            tweetable_times = []
            for node in tweetable_nodes:
                tweetable_times.extend([node.created_at, node.updated_at, node.published_at])

            new_nodes = []

            # this will rate limit
            # if something goes wrong, it will not duplicate tweets
            tweet_nodes(tweetable_nodes, start_time, dry_run=False)

            start_time = datetime.now(TIMEZONE)

        time.sleep(3600)


if __name__ == "__main__":
    main()
