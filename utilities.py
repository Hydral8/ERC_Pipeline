import asyncio
from collections import namedtuple

import aiohttp
import itertools
import math
import os
import os.path as osp
import sys
import traceback
from typing import Tuple, Dict, List, Callable, Any, Union

import xlsxwriter
from ete3 import PhyloTree
from ete3.parser.newick import NewickError
import networkx as nx
import csv
from statsmodels.stats.multitest import multipletests
import numpy as np
import pandas as pd


TaxaInfo = namedtuple('TaxaInfo', ['txid', 'name', 'species', 'genus', 'family', 'order'])
Letters = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"]

color_dict = {
    "yellow": "#cae32b",
    "orange": "#edb86d"
}

def safe_phylo_read(filename) -> PhyloTree:
    if isinstance(filename, PhyloTree):
        return filename
    try:
        return PhyloTree(filename, format=3)
    except:
        try:
            return PhyloTree(filename)
        except:
            try:
                return PhyloTree(filename, format=1)
            except:
                try:
                    return PhyloTree(filename, format=5)
                except NewickError as e:
                    print(f"Are you sure tree {filename} exists?", file=sys.stderr, flush=True)
                    raise e


def mammal_taxa_info(name_as_key: bool = False) -> Dict[str, TaxaInfo]:
    info = dict()
    with open(osp.join(_self_path(), 'data', 'txid2name.tsv'), 'r') as f:
        first = True
        for line in f:
            if first:
                first = False
                continue
            split = line.strip().split("\t")
            info[split[1] if name_as_key else split[0]] = TaxaInfo(*split)
    return info


def safe_delete(path: str):
    if osp.exists(path):
        try:
            os.remove(path)
        except: pass


def safe_mkdir(dir: str):
    if not osp.exists(dir):
        try:
            os.makedirs(dir)
        except: pass


def wait(coro):  # Wait on coroutines more simply
    return asyncio.get_event_loop().run_until_complete(coro)


# noinspection PyPep8Naming
def translate(G: nx.Graph, l2n: Dict[str, str]) -> nx.Graph:
    G_c = nx.Graph() if not G.is_directed() else nx.DiGraph()

    for u, v in G.edges:
        G_c.add_edge(l2n.get(u, u), l2n.get(v, v), **G.get_edge_data(u, v), default=dict())
    return G_c


def override_sys_out(tag: str = None):
    import sys
    from datetime import datetime as dt

    def new_write(iostream):

        orig = iostream.write

        class TimeStamper:

            nl = True

            def write(self, s):
                """Write function overloaded."""
                if s == '\n':
                    orig(s)
                    self.nl = True
                elif self.nl:
                    orig('[%s]%s: %s' % (str(dt.now()), '' if not tag else f"[{tag}]", s))
                    self.nl = False
                else:
                    orig(s)

        stamper = TimeStamper()

        iostream.write = stamper.write

    new_write(sys.stdout)
    new_write(sys.stderr)


async def async_call(cmd, cwd=None, **kwargs):
    print("> " + cmd, flush=True)
    if "shell" in kwargs.keys():
        del kwargs["shell"]
    kwargs['cwd'] = cwd if cwd else os.getcwd()

    process = await asyncio.create_subprocess_shell(cmd, **kwargs)
    await process.wait()


def chunks(l, n):
    # For item i in a range that is a length of l,
    for i in range(0, len(l), n):
        # Create an index range for l of n items:
        yield l[i:i+n]


def try_hook_uvloop():  # Can improve speed on unix systems
    try:
        import uvloop
        uvloop.install()
    except: pass


def _self_path():
    path = osp.dirname(__file__)
    if not path:
        path = '.'
    return path


def _add_edge(net: nx.Graph, a, b):
    if not net.has_edge(a, b):
        net.add_edge(a, b)


def rho_sorted_neighbors(net: nx.Graph, node: str) -> List[str]:
    return list(sorted(net.neighbors(node), key=lambda n: net[n][node]['rho'], reverse=True))

def bonferroni_sorted_neighbors(net: nx.Graph, node: str, protein_list: List[str]) -> List[str]:
    neighbors = list(net.neighbors(node))
    proteins_of_interest = neighbors
    if protein_list is not None:
        proteins_of_interest = [prot for prot in protein_list if prot in set(neighbors)]
        
    n_neighbors = len(neighbors)
    bonferroni_p = [net[node][n]['p'] * float(n_neighbors) for n in proteins_of_interest]
    bonferroni_p_dict = {n: bonferroni_p[i] for i,n in enumerate(proteins_of_interest)}
    return list(sorted(proteins_of_interest, key=lambda n: bonferroni_p_dict[n]))

def get_ranks(sorted:list):
    return {node: i + 1 for i,node in enumerate(sorted)}


def p_sorted_neighbors(net: nx.Graph, node: str) -> List[str]:
    return list(sorted(net.neighbors(node), key=lambda n: net[n][node]['p']))

def sorted_dict(net: nx.Graph, node: str, type: str) -> dict:
    neighbours = net.neighbors(node)
    # test if self works too (i.e. net[n][n][p] and change code accordingly if it does)
    v_dict = {n: net[n][node][type] for n in neighbours}
    
    return dict(sorted(v_dict.items(), key = lambda k: k[1]))

def rho_matrix(net: nx.Graph, filename: str, l2n: dict = None, cutoff=0.9, operator=">", protein_list=None, no_csv=False):
    nodes = list(sorted(net.nodes), key = lambda n: l2n.get(n, n))
    headers = [l2n[prot] for prot in nodes] if l2n is not None else nodes
    # excel sheet
    book = xlsxwriter.Workbook(f"{filename}.xlsx")
    rho_sheet = book.add_worksheet("rho_matrix")
    #formatting
    bold = make_bold_formatting(book)
    rho_format = make_rho_formatting(book)
    highlight = make_highlight_formatting(book)
    
    headers.insert(0, "ROI")
    
    rho_sheet.write_row(0, 1, [l2n[prot] for prot in nodes] if l2n is not None else nodes, cell_format=bold)
    rho_sheet.write_column(0, 0, headers, cell_format=bold)
    
    data = []
    
    for i, id1 in enumerate(nodes):
        row_data = [id1]
        for j, id2 in enumerate(nodes):
            if id1 == id2:
                row_data.append(1.0)
            else:
                # print(id1 + " " + id2)
                try:
                    row_data.append(net[id1][id2]['rho'])
                except Exception as e:
                    print("Failed Exporting erc data between : " + id1 + " " + id2)
        rho_sheet.write_row(i + 1, 1, row_data[1:])
        data.append(row_data)
        
    col_name = "B"
    max_row = max_col = len(nodes)
    
    
    rho_sheet.conditional_format(1,1, max_row, max_col, {
        'type': 'formula',
        'criteria': f"=AND({col_name}2 {operator} {cutoff}, NOT(ISBLANK({col_name}2)))",
        'format': highlight
    })
    
        
    book.close()
        
    if not no_csv:
        with open(f'{filename}.csv', 'w', encoding='UTF8') as f:
            writer = csv.writer(f)
            
            writer.writerow(headers)
            
            for row_data in data:
                writer.writerow(row_data)
    

def p_matrix(net: nx.Graph, filename: str, l2n: dict = None, correction="bonferroni", cutoff=0.01, operator="<", protein_list=None, no_csv=False):
    nodes = list(net.nodes)
    headers = [l2n[prot] for prot in list(net.nodes)] if l2n is not None else list(net.nodes)
    # excel sheet
    book = xlsxwriter.Workbook(f"{filename}.xlsx")
    p_sheet = book.add_worksheet(f"{correction}_p_matrix")
    #formatting
    bold = make_bold_formatting(book)
    p_format = make_p_formatting(book)
    highlight = make_highlight_formatting(book)
    
    headers.insert(0, "ROI")
    
    p_sheet.write_row(0, 1, [l2n[prot] for prot in nodes] if l2n is not None else nodes, cell_format=bold)
    p_sheet.write_column(0, 0, headers, cell_format=bold)
    
    data = []
    
    n_neighbors = len(nodes) - 1
    
    for i, id1 in enumerate(nodes):
        row_data = [id1]
        for j, id2 in enumerate(nodes):
            if id1 == id2:
                row_data.append(1.0)
            else:
                # print(id1 + " " + id2)
                try:
                    if correction == "bonferroni":
                        row_data.append(net[id1][id2]['p'] * n_neighbors)
                    else:
                        row_data.append(net[id1][id2]['p'])
                except Exception as e:
                    print("Failed Exporting erc data between : " + id1 + " " + id2)
        p_sheet.write_row(i + 1, 1, row_data[1:])
        data.append(row_data)
        
    col_name = "B"
    max_row = max_col = len(nodes)
    
    
    p_sheet.conditional_format(1,1, max_row, max_col, {
        'type': 'formula',
        'criteria': f"=AND({col_name}2 {operator} {cutoff}, NOT(ISBLANK({col_name}2)))",
        'format': highlight
    })
    
    book.close()
    
    if not no_csv:
        with open(f'{filename}.csv', 'w', encoding='UTF8') as f:
            writer = csv.writer(f)
            
            writer.writerow(headers)
            
            for row_data in data:
                writer.writerow(row_data)
    
    

def rank_matrix(net: nx.Graph, filename: str, l2n: dict = None):
    nodes = list(net.nodes)
    headers = [l2n[prot] for prot in list(net.nodes)] if l2n is not None else list(net.nodes)
    # excel sheet
    book = xlsxwriter.Workbook(f"{filename}.xlsx")
    rank_sheet = book.add_worksheet("ranks")
    #formatting
    bold = make_bold_formatting(book)
    rank_format = make_rank_formatting(book)
    
    headers.insert(0, "POI")
    
    rank_sheet.write_row(0, 1, [l2n[prot] for prot in nodes] if l2n is not None else nodes, cell_format=bold)
    rank_sheet.write_column(0, 0, headers, cell_format=bold)
    
    row_idx = 0
    
    with open(f'{filename}.csv', 'w', encoding='UTF8') as f:
        writer = csv.writer(f)
        
        writer.writerow(headers)
        
        for node in nodes:
            row_idx += 1
            ranks = get_ranks(bonferroni_sorted_neighbors(net, node))
            rank_list = [l2n[node] if l2n is not None else node]
            for subnode in nodes:
                if subnode in ranks:
                    rank_list.append(ranks[subnode])
                elif subnode == node:
                    rank_list.append(1)
                else:
                    rank_list.append(0)
            writer.writerow(rank_list)
            rank_sheet.write_row(row_idx,1, rank_list[1:], cell_format=rank_format)
            
    book.close()

def get_col_letter(idx):
    if idx > 25:
        quotient = (idx // 26) - 1
        remainder = idx % 26
        return Letters[quotient] + Letters[remainder]
    else:
        return Letters[idx]
    
    
    
    
        
def rank_protein(net: nx.Graph, prots: List[str], filename: str, l2n: dict = None, datatype='bonferroni_p', cutoff=0.01, operator="<", protein_list=None, use_protein_list_only = True, no_csv=False):
    
    """
    TODO: convert get_ranks to rank dict so that we can just reference that rank dict when we want to get the rank of a partner in a protein and protein in a partner. 
    Single store of ranks rather than generating ranks multiple times for each protein through a loop
    """
    
    dirpath = os.path.dirname(filename)
    
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)
        
    #check if protein list exists -> if it does if its a path open it and make into list
    if isinstance(protein_list, str):
        with open(protein_list, "r") as f:
            protein_list = list(map(lambda k: k.strip(), f.read().splitlines()))

    book = xlsxwriter.Workbook(f"{filename}.xlsx")
    rank_sheet = book.add_worksheet("ranks")
    #formatting
    bold = make_bold_formatting(book)
    rank_format = make_rank_formatting(book)
    highlight_blue = make_highlight_formatting(book)
    highlight_orange = make_highlight_formatting(book, color="orange")
    highlight_yellow = make_highlight_formatting(book, color="yellow")
    
    shift = 3
    
    if protein_list is not None:
        rank_sheet.write_column(0, 0, ["protein_list", *list(map(lambda k: l2n[k], protein_list))])
    else:
        rank_sheet.write_column(0, 0, ["protein_list", "ALL"])
    
    for prot in prots:
        
        if use_protein_list_only:
            sorted = bonferroni_sorted_neighbors(net, prot, protein_list)
        else:
            sorted = bonferroni_sorted_neighbors(net, prot, None)
        
        rank_in_partners = []
        
        if protein_list and (prot not in protein_list):
            protein_list.append(prot)
        
        for sorted_prot in sorted:
            if use_protein_list_only:
                rank_in_partners.append(get_ranks(bonferroni_sorted_neighbors(net, sorted_prot, protein_list))[prot])
            else:
                rank_in_partners.append(get_ranks(bonferroni_sorted_neighbors(net, sorted_prot, None))[prot])
            
        n_neighbors = len(list(net.neighbors(prot)))
        bonferroni_p = [net[prot][n]['p'] * float(n_neighbors) for n in sorted]
        bonferroni_p_dict = {n: bonferroni_p[i] for i,n in enumerate(sorted)}
            
        # fdr_corrected = multipletests([net[n][n2]['p'] for n2 in node2sorted[n]], method="fdr_bh")[1]
        headers = ["prot", f"rank in {l2n[prot] if l2n is not None else prot}", f"{l2n[prot] if l2n is not None else prot} rank in partner", "rho", "p", "bonferroni_p"]
        proteins = np.array([l2n[sorted_prot] for sorted_prot in sorted] if l2n is not None else sorted)
        ranks = np.array([i + 1 for i,x in enumerate(sorted)])
        rhos = np.array([net[prot][n]['rho'] for n in sorted])
        ps = np.array([net[prot][n]['p'] for n in sorted])
        bonferroni_ps = np.array([bonferroni_p_dict[n] for n in sorted])
        
        rank_sheet.write_row(0, shift + 0, headers, cell_format=bold)
        rank_sheet.write_column(1, shift + 0, proteins, cell_format=bold)
        rank_sheet.write_column(1, shift + 1, ranks, cell_format=bold)
        rank_sheet.write_column(1, shift + 2, rank_in_partners, cell_format=bold)
        rank_sheet.write_column(1, shift + 3, rhos, cell_format=bold)
        rank_sheet.write_column(1, shift + 4, ps, cell_format=bold)
        rank_sheet.write_column(1, shift + 5, bonferroni_ps, cell_format=bold)
        
        # select col based on datatype
        cols = {
            "prot": [shift, get_col_letter(shift)],
            "rank": [shift + 1, get_col_letter(shift + 1)],
            "rank_partner": [shift + 2, get_col_letter(shift + 2)],
            "rho": [shift + 3, get_col_letter(shift + 3)],
            "p": [shift + 4, get_col_letter(shift + 4)],
            "bonferroni_p": [shift + 5, get_col_letter(shift + 5)]
        }
        
        max_row = len(proteins)
        col_num = cols[datatype][0]
        col_name = cols[datatype][1]
        
        shift += 7
        
        rank_sheet.conditional_format(1,col_num, max_row, col_num, {
            'type': 'formula',
            'criteria': f"=AND({col_name}2 {operator} {cutoff}, NOT(ISBLANK({col_name}2)))",
            'format': highlight_blue
        })
        
        
        rank_col_num = cols["rank"][0]
        rank_col_name = cols["rank"][1]
        rank_partner_col_num = cols["rank_partner"][0]
        rank_partner_col_name = cols["rank_partner"][1]
        rank_sheet.conditional_format(1,rank_col_num, max_row, rank_partner_col_num, {
            'type': 'formula',
            'criteria': f"=AND(AND({rank_col_name}2 <= 30, {rank_partner_col_name}2 <= 30), NOT(ISBLANK({rank_col_name}2)))",
            'format': highlight_orange
        })
        
        if (protein_list is not None) and (use_protein_list_only is False):
        
            col_num = cols["prot"][0]
            col_name = cols["prot"][1]
            prot_list_length = len(protein_list)
            
            rank_sheet.conditional_format(1,col_num, max_row, col_num, {
                'type': 'formula',
                'criteria': f"=AND(COUNTIF($A$2:$A${prot_list_length}, {col_name}2) > 0, NOT(ISBLANK({col_name}2)))",
                'format': highlight_yellow
            })
        
        
        
    
        if not no_csv:
            name = l2n[prot]
            empty = np.empty(len(bonferroni_ps)).fill("")
            data = np.stack([prots, ranks, rank_in_partners, rhos, ps, bonferroni_ps, empty], -1)
            
            with open(f"{filename}_{name}.csv", "w") as f:
                f.write(",".join(headers) + "\n")
                for vals in data:
                    f.write(",".join(vals) + "\n")
                    
    book.close()
        
# TODO: ADD COMMENTS
def compare(filename1: str, filename2: str, outputfile: str, l2n: dict = None):
    book = xlsxwriter.Workbook(f"{outputfile}.xlsx")
    data_sheet = book.add_worksheet("data")
    
    headers_base = None
    data_vals_first = []
    name_dict = dict()
    padding = 3
    bold = make_bold_formatting(book)
    rank_format = make_rank_formatting(book)
    
    k = 0
    with open(filename1, "r") as f:
        it = iter(f)
        headers_base = next(it).strip().split(",")
        name = headers_base[1].split()
        name[-1] = l2n[name[-1]] if l2n else name[-1]
        headers_base[1] = (" ").join(name)
        name = headers_base[2].split()
        name[0] = l2n[name[0]] if l2n else name[0]
        headers_base[2] = (" ").join(name)
        # configure naming with l2n -> switch to known protein names rather than number values
        
        for l in it:
            vals = l.strip().split(",")
            vals[0] = l2n[vals[0]] if l2n else vals[0]
            name_dict[vals[0]] = None

            k += 1
            # add padding
            for i in range(padding):
                vals.append("")
            data_vals_first.append(vals)
            
    data_vals_second = []
    with open(filename2, "r") as f:
        it = iter(f)
        next(it) # skip header
        for l in it:
            vals = l.strip().split(",")
            vals[0] = l2n[vals[0]] if l2n else vals[0]
            name_dict[vals[0]] = vals
            vals_padded = vals.copy()
            for i in range(padding):
                vals_padded.append("")
            data_vals_second.append(vals_padded)
    
    sorted_vals_second = []
    for (key, vals) in name_dict.items():
        if vals is not None:
            for i in range(padding):
                vals.append("")
            sorted_vals_second.append(vals)
        
    data = np.concatenate([data_vals_first, data_vals_second, data_vals_first, sorted_vals_second], axis=-1)
    
    types = []
    types_base = [('prot', '<U16'), ('rank', 'i'), ('rank_in_partners','i'), ('rhos', 'f'), ('ps', 'f'), ('bonferroni_ps', 'f')]
    for i in range(4):
        for k in types_base:
            types.append((f"{k[0]}_{i}", k[1]))
        for j in range(padding):
            types.append((f'space{j}_{i}', '<U1'))

    data = np.array(list(map(tuple, data)), dtype=types)
    
    headers = []
    for i in range(4):
        for k in headers_base:
            headers.append(k)
        for j in range(padding):
            headers.append("")
    data_sheet.write_row(0, 1, headers, cell_format=bold)
    idx = 1
    for vals in data:
        data_sheet.write_row(idx, 1, vals)
        idx += 1
    
    book.close()
    
    

def rrn_nets(net: nx.Graph, *odbs: str, types=["steps"]) -> Tuple[nx.DiGraph, nx.DiGraph, nx.DiGraph]:
    step1 = nx.DiGraph()
    added_set = set()
    base_10_added_set = set()
    types = set(types)
    
    final_return = {}

    # Reciprank 20 first
    for prot in odbs:
        for i, neighbor in enumerate(rho_sorted_neighbors(net, prot)[:20]):
            prot_pos = rho_sorted_neighbors(net, neighbor).index(prot)
            if prot_pos < 20:
                added_set.add(neighbor)
                _add_edge(step1, prot, neighbor)
                _add_edge(step1, neighbor, prot)
                
    if "base" in types:
        final_return["base"] = step1
    # reciprank RR10
    if "base10" in types:
        base_10 = nx.DiGraph()
        for prot in odbs:
            for i, neighbor in enumerate(rho_sorted_neighbors(net, prot)[:10]):
                prot_pos = rho_sorted_neighbors(net, neighbor).index(prot)
                if prot_pos < 20:
                    base_10_added_set.add(neighbor)
                    _add_edge(base_10, prot, neighbor)
                    _add_edge(base_10, neighbor, prot)
        final_return["base10"] = base_10
    if "steps" in types:
        # Reciprank 20 of added
        step2 = step1.copy()
        for added in added_set:
            for i, neighbor in enumerate(rho_sorted_neighbors(net, added)[:20]):
                prot_pos = rho_sorted_neighbors(net, neighbor).index(added)
                if prot_pos < 20:
                    _add_edge(step2, added, neighbor)
                    _add_edge(step2, neighbor, added)

        # Fill connections <=20 between current node set
        step3 = step2.copy()
        for (n1, n2) in itertools.permutations(step2.nodes, 2):
            if n1 == n2 or step2.has_edge(n1, n2):
                continue
            prot_pos = rho_sorted_neighbors(net, n1).index(n2)
            if prot_pos < 20:
                _add_edge(step3, n1, n2)
        
        final_return["step2"] = step2
        final_return["step3"] = step3
        
    if "closed":
        closed_nx = nx.DiGraph()
        for prot in odbs:
            for i, neighbor in enumerate(rho_sorted_neighbors(net, prot)[:20]):
                if neighbor in odbs:
                    _add_edge(closed_nx, prot, neighbor)
        final_return["closed"] = closed_nx
        
    if "closed10":
        closed_nx = nx.DiGraph()
        for prot in odbs:
            for i, neighbor in enumerate(rho_sorted_neighbors(net, prot)[:10]):
                if neighbor in odbs:
                    _add_edge(closed_nx, prot, neighbor)
        final_return["closed10"] = closed_nx
        
    # all of step1 base network + any connections between partners of central protein
    # Add unidirectional 20 (and thereby reciprocal 20's) of each protein within protein set. Only adding edges not new proteins
    if "closed_base" in types:   
        closed_nx = step1.copy()
        for added in added_set:
            for i, neighbor in enumerate(rho_sorted_neighbors(net, added)[:20]):
                if neighbor in added_set:
                    _add_edge(closed_nx, added, neighbor)
        final_return["closed_base"] = closed_nx
    
    if "closed_base10" in types:
        closed_nx = base_10.copy()
        for added in base_10_added_set:
            for i, neighbor in enumerate(rho_sorted_neighbors(net, added)[:10]):
                if neighbor in base_10_added_set:
                    _add_edge(closed_nx, added, neighbor)
        final_return["closed_base10"] = closed_nx
        
    # if original central protein (not partners) in top 20 unidirecitonal
    if "prot_in_20" in types:
        prot_in_20 = nx.DiGraph()
        for prot in odbs:
            neighbors = rho_sorted_neighbors(net, prot)
            for node in neighbors:
                if rho_sorted_neighbors(net, node).index(prot) < 20:
                    _add_edge(prot_in_20, node, prot)
        final_return["prot_in_20"] = prot_in_20
        
    if "prot_in_10" in types:
        prot_in_10 = nx.DiGraph()
        for prot in odbs:
            neighbors = rho_sorted_neighbors(net, prot)
            for node in neighbors:
                if rho_sorted_neighbors(net, node).index(prot) < 10:
                    _add_edge(prot_in_10, node, prot)
        final_return["prot_in_10"] = prot_in_10
        
    return final_return


def make_rho_formatting(workbook: xlsxwriter.Workbook):
    return workbook.add_format({
        "num_format": "0.000"
    })


def make_p_formatting(workbook: xlsxwriter.Workbook):
    return workbook.add_format({
        "num_format": "0.00E+00"
    })
    
def make_rank_formatting(workbook: xlsxwriter.Workbook):
    return workbook.add_format({
        "num_format": "0"
    })


def make_bold_formatting(workbook: xlsxwriter.Workbook):
    return workbook.add_format({
        'bold': True
    })
    
def make_highlight_formatting(workbook: xlsxwriter.Workbook, color="#13cdd4"):
    if "#" not in color:
        color = color_dict[color]
    
    return workbook.add_format({
        'bg_color': color
    })

def color_from_custom_map(pct: float, cmap: Dict[float, Tuple[int, int, int]]) -> str:
    # https://stackoverflow.com/a/7128796/5179044
    maxColor = 0
    lastMaxColor = 0
    for col in sorted(cmap.keys()):
        lastMaxColor = maxColor
        maxColor = col
        if pct < col:
            break

    _range = maxColor - lastMaxColor
    range_pct = (pct - lastMaxColor) / _range

    pct_lower = 1 - range_pct
    pct_upper = range_pct

    color = (
        math.floor(cmap[lastMaxColor][0] * pct_lower + cmap[maxColor][0] * pct_upper),
        math.floor(cmap[lastMaxColor][1] * pct_lower + cmap[maxColor][1] * pct_upper),
        math.floor(cmap[lastMaxColor][2] * pct_lower + cmap[maxColor][2] * pct_upper),
    )

    hex = "#{:02x}{:02x}{:02x}".format(color[0], color[1], color[2])

    return hex


def graphviz_network_plot(net: nx.Graph, output: str, highlight: Dict[str, str] = dict(), circo: bool = False,
                          highlight_by_font: bool = False):
    # output: a png file
    # Requires pygraphviz

    cmap = {0.0: (255, 255, 255),
            0.5: (204, 217, 255),
            1.0: (20, 139, 255)}

    max_connectivity = max([len(list(net.neighbors(n))) for n in net.nodes])

    A = nx.nx_agraph.to_agraph(net)
    A.graph_attr.update(strict=False, overlap=False, splines='true')
    A.node_attr['style'] = 'filled'
    if circo:
        A.layout(prog="circo")
    else:
        A.layout()
    for n in net.nodes:
        n2 = A.get_node(n)
        connectivity = len(list(net.neighbors(n)))
        n2.attr["fillcolor"] = color_from_custom_map(connectivity / max_connectivity, cmap) if highlight_by_font or (not highlight_by_font and n not in highlight) else highlight[n]
        if highlight_by_font and n in highlight:
            n2.attr['fontcolor'] = highlight[n]
    A.draw(output, args='-Gsize=10,10 -Gdpi=300')
    A.draw(output.replace(".png", ".svg"), format='svg')


def decorate(fun: Callable[..., Any]) -> type:
    class decorator:

        def __init__(self):
            pass

        def __call__(self, f):
            def wrapped(*args, **kwargs):
                try:
                    o = fun(*args, **kwargs)
                    if o is not None:
                        return o
                except:
                    pass
                return f(*args, **kwargs)

            return wrapped

    return decorator


RECURSION_LIMIT = 5


class Limiter:

    def __init__(self, delay: int = 1):
        self._time = -1
        self._limiter_lock = asyncio.Lock()
        self.request_delay = delay
        self.changed = False

    async def __aenter__(self):
        await self._limiter_lock.acquire()

        to_wait = self._time_to_wait()

        await asyncio.sleep(to_wait)

    async def __aexit__(self, *args, **kwargs):
        self._update_timer()

        self._limiter_lock.release()

    def _time_to_wait(self) -> int:
        request_time = self._time

        if request_time == -1:
            return 0

        now = asyncio.get_event_loop().time()

        to_wait = request_time + self.request_delay - now
        to_wait = max(0, to_wait)

        return to_wait

    def _update_timer(self):
        now = asyncio.get_event_loop().time()
        self._time = now


_new_session_lock = asyncio.Lock()
_limiters = dict()

_sessions = dict()


async def remote_call(url: str, is_json: bool, *, recursion_count: int = 0, **kwargs: Any) -> Union[Dict, str]:
    orig_url = url
    if recursion_count >= RECURSION_LIMIT:
        raise Exception(f"Too much recursion! (url={url})")

    if kwargs is not None and len(kwargs) > 0:
        first = True
        for k, v in kwargs.items():
            if v is None:
                continue

            if first:
                url = url + '?'
                first = False
            else:
                url = url + '&'

            if isinstance(v, set) or isinstance(v, list):
                v = ",".join(v)

            url = url + f'{k}={v}'

    base_url = url.split('/')[0]

    await _new_session_lock.acquire()

    if base_url not in _sessions:
        _sessions[base_url] = aiohttp.ClientSession()
        _limiters[base_url] = Limiter()

    session = _sessions[base_url]
    limiter = _limiters[base_url]

    _new_session_lock.release()

    async with limiter:
        try:
            async with session.get(url) as resp:
                resp.raise_for_status()
                if not limiter.changed and resp.headers.get('x-rate-limit-limit', None) is not None:  # Prepare ratelimiting if applicable
                    req_limit = float(resp.headers.get('x-rate-limit-limit'))  # We are assuming a 60 second window
                    limit_time = float(resp.headers.get('x-rate-limit-reset', "1"))  # We are assuming that the first request will give us the window time
                    limiter.request_delay = limit_time / req_limit
                    limiter.changed = True
                if resp.status != 503:
                    if resp.status == 429:  # Rate limit
                        if int(resp.headers.get('x-rate-limit-remaining', 0)) == 0:
                            await asyncio.sleep(float(resp.headers.get('x-rate-limit-reset', "1")))
                    elif is_json:
                        return await resp.json(content_type=None)
                    else:
                        return await resp.text()
        except Exception as e:
            print(f"Exception for url {url} caught!", flush=True)
            traceback.print_exc()
            await asyncio.sleep(1)  # Sleep a little

    return await remote_call(orig_url, is_json, recursion_count=recursion_count + 1, **kwargs)


async def close_web_sessions():
    for session in _sessions.values():
        await session.close()
    _sessions.clear()
