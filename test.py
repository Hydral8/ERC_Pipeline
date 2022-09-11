from pipeline import *
from utilities import *
import asyncio
import numpy
from statsmodels.stats.multitest import multipletests

def main():
    register_previous_run("../../../research/werren_lab/new_ace2_with_old_dataset")
    erc_network = generate_network()
    l2n = make_l2n(True)
    ace2 = "59099at40674"
    step1, step2, step3 = rrn_nets(erc_network, ace2)
    step2 = translate(step2, l2n)
    
    graphviz_network_plot(step3, "step3.png")
    graphviz_network_plot(step2, "step2.png")
    # rank_matrix(erc_network, "./Analysis_Results/old_ace2_new_dataset/ranks", l2n)
    # rank_protein(erc_network, "59099at40674", "./Analysis_Results/old_ace2_new_dataset/rank-single_no_l2n")
    # rank_protein(erc_network, "59099at40674", "./Analysis_Results/new_ace2_old_dataset/rank-single", l2n)


def main_2():
    p_vals = []
    data_vals = []
    bonferroni_p = []
    with open("ACE2_to_and_from_lists.csv", "r") as f:
        it = iter(f)
        headers = next(it).strip().split(",")
        for l in it:
            vals = np.array(l.strip().split(","))
            data_vals.append(vals)
            p_vals.append(float(vals[-2])) 
        bonferroni_p = multipletests(p_vals, method="bonferroni")[1]
        
    data = np.array(data_vals)
    data[:, -1] = bonferroni_p
    
    with open("ACE2_to_and_from_lists.csv", "w") as f:
        f.write(",".join(headers) + "\n")
        for vals in data:
            f.write(",".join(vals) + "\n")
            
def main_3():
    l2n = make_l2n(True)
    compare("./analysis_results/old_ace2_new_dataset/rank-single_no_l2n.csv", "./analysis_results/new_ace2_new_dataset/rank-single_no_l2n.csv", "./analysis_results/compare/old_vs_new_ace2_new_dataset", l2n=l2n)
            
def main_4():
    register_previous_run("../../../research/werren_lab/results/RTL9/10kmer_sliding_window_results")
    erc_network = generate_network()
    p_matrix(erc_network, "./Analysis_Results/RTl9_10kmer_sliding_window/p_matrix")

def main_5():
    register_previous_run("../../../research/werren_lab/new_ace2_with_old_dataset")
    erc_network = generate_network()
    l2n = make_l2n(True)
    rtl9 = "35692at40674"
    # step1, step2, step3 = rrn_nets(erc_network, rtl9)
    # closed_net, = rrn_nets(erc_network, rtl9, types=["closed"])
    prot_in_20_net, = rrn_nets(erc_network, rtl9, types=["prot_in_20"])
    
    # grap top 50 proteins of rtl9
    neighbours = p_sorted_neighbors(erc_network, rtl9)
    neighbours_100 = neighbours[:100]
    neighbours_200 = neighbours[:200]
    
    # rank_rtl9_in_20 = []
    # for prot in neighbours:
    #     if p_sorted_neighbors(erc_network, prot).index(rtl9) < 20:
    #         rank_rtl9_in_20.append(prot)
            
    # enrich_network(rank_rtl9_in_20, erc_network, l2n, "./analysis_results/rtl9/enrichments/rank_rtl9_in_20_enrichments")
    # enrich_network(neighbours_100, erc_network, l2n, "./analysis_results/rtl9/enrichments/top_100_enrichments")
    # enrich_network(neighbours_200, erc_network, l2n, "./analysis_results/rtl9/enrichments/top_200_enrichments")
    # closed_net = translate(closed_net, l2n)
    prot_in_20_net = translate(prot_in_20_net, l2n)
    
    graphviz_network_plot(prot_in_20_net, "./analysis_results/rtl9/prot_in_20_net.png")
    
def main_6():
    register_previous_run("../../../research/werren_lab/analysis/new_ace2_with_old_dataset")
    register_previous_run("../../../research/werren_lab/analysis/fu_plus_1953")
    erc_network = generate_network()
    l2n = make_l2n(True)
    
    #prots = ["35956at40674","66939at40674", "71150at40674", "28286at40674", "854at40674", "154212at40674"] #elp2, elp3, TRMT1, nat10, kmt2d, METTL1
    with open("./protein_lists/fu.txt", "r") as f:
        prots = list(map(lambda k: k.strip(), f.read().splitlines()))
    prots.append("854at40674")
    
    # for prot in prots:
    rank_protein(erc_network, prots, f"../../../research/werren_lab/analysis/fu_plus_1953/fu_proteins_combined/rank-single_with_l2n", l2n, protein_list="./protein_lists/fu.txt", use_protein_list_only = False, no_csv=True)
    # rank_protein(erc_network, prots, f"../../../research/werren_lab/analysis/fu_plus_1953/fu_proteins_combined/rank-single_with_l2n_only_fu", l2n, protein_list="./protein_lists/fu.txt", no_csv=True)
    # rank_protein(erc_network, "59099at40674", "./Analysis_Results/new_ace2_old_dataset/rank-single", l2n)
    
def main_7():
    register_previous_run("../../../research/werren_lab/analysis/new_ace2_with_old_dataset")
    register_previous_run("../../../research/werren_lab/analysis/fu_plus_1953")
    erc_network = generate_network()
    l2n = make_l2n(True)
    
    p_matrix(erc_network, "../../../research/werren_lab/analysis/fu_plus_1953/bonferroni_corrected_p_matrix", l2n, no_csv=True)
    
def main_8():
    register_previous_run("../../../research/werren_lab/analysis/new_ace2_with_old_dataset")
    register_previous_run("../../../research/werren_lab/analysis/fu_plus_1953")
    erc_network = generate_network()
    l2n = make_l2n(True)
    
    with open("./protein_lists/fu.txt", "r") as f:
        prots = list(map(lambda k: k.strip(), f.read().splitlines()))
    # prots.append("854at40674") KMT2D
    
    data_dict = rrn_nets(erc_network, *prots, types=["base10", "closed10", "closed", "closed_base10", "prot_in_10"])
    
    base10 = data_dict["base10"]
    closed = data_dict["closed"]
    closed10 = data_dict["closed10"]
    closed_base10 = data_dict["closed_base10"]
    prot_in_10 = data_dict["prot_in_10"]
    
    base10 = translate(base10, l2n)
    closed = translate(closed, l2n)
    closed10 = translate(closed10, l2n)
    closed_base10 = translate(closed_base10, l2n)
    prot_in_10 = translate(prot_in_10, l2n)
    
    
    graphviz_network_plot(base10, "./analysis_results/dragony_fu/base10.png")
    graphviz_network_plot(closed, "./analysis_results/dragony_fu/closed.png")
    graphviz_network_plot(closed10, "./analysis_results/dragony_fu/closed10.png")
    graphviz_network_plot(closed_base10, "./analysis_results/dragony_fu/closed_base10.png")
    graphviz_network_plot(prot_in_10, "./analysis_results/dragony_fu/prot_in_10.png")
    
     
    
if __name__ == "__main__":
    # main()
    # main_2()  
    # # main_3()    
    # main() 
    # main_4()
    # main_5()
    main_6()
    # main_7()
    # main_8()
