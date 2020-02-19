#!/usr/bin/env python3

from Bio import Entrez
import xmltodict
from box import Box
from datetime import datetime


def get_nodes() -> Box:
    """
    Returns a list of nodes that have been created/updated/published after a
    specified date
    """

    start_time = datetime.now().strftime("%Y/%m/%d")

    # I haven't been able to figure out how to query by hour so we have to
    # parse the times and filter manually. No biggie.
    term = f'"{start_time}"[EDAT] : "3000"[EDAT]'

    node_list = Box(
        xmltodict.parse(Entrez.esearch("taxonomy", term, email="austin@onecodex.com").read())
    )

    tax_ids = node_list.eSearchResult.IdList.Id[:2]

    nodes = Box(
        xmltodict.parse(
            Entrez.efetch(db="taxonomy", id=",".join(tax_ids), email="austin@onecodex.com").read()
        )
    )

    return nodes.TaxaSet.Taxon


nodes = get_nodes()

for node in nodes:
    print(node.ScientificName)
    print([node.CreateDate, node.UpdateDate, node.PubDate])
