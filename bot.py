#!/usr/bin/env python3

from Bio import Entrez
import xmltodict
from box import Box
from datetime import datetime
from dateutil.parser import parse as parse_date
from pprint import pprint
from typing import List
import os
import json


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


def format_tweet_for_node(node) -> str:
    if "Bacteria" in node.lineage:
        emoji = "🦠"
    elif "Fungi" in node.lineage:
        emoji = "🍄"
    else:
        emoji = ""

    return "\n".join([f"{node.name} ({node.rank}) {emoji}".strip(), f"lineage: {node.lineage}"])


def main():

    if os.path.exists("last-time.json"):
        with open("last-time.json") as handle:
            last_time = parse_date(json.load(handle)["last-time"])
    else:
        last_time = datetime.now()

    with open("last-time.json", "w") as handle:
        json.dump({"last-time": datetime.now()}, handle, default=str)

    print(f"looking up taxa that were created since {last_time}")

    nodes = get_nodes(start_time=last_time)

    print(f"got {len(nodes)} new nodes")

    new_nodes = [n for n in nodes if (n.created_at >= start_time) or (n.published_at >= start_time)]
    new_nodes = nodes

    for node in new_nodes:
        print(format_tweet_for_node(node))
        print()


if __name__ == "__main__":
    main()
