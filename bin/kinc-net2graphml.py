#!/usr/bin/env python3
"""
Converts a KINC network into an GraphML file

"""
import os
import argparse
import numpy as np
import pandas as pd
import igraph as ig
from fa2 import ForceAtlas2


def load_network(file_path):
    """
    Imports the KINC-generated network file (either full or Tidy versions).

    file_path : The path to the network file.

    return : A pandas dataframe containing the network.
    """

    net = pd.read_csv(file_path, sep="\t")

    # Make sure the file has the required columns
    columns = net.columns
    if ('Source' not in columns) | ('Target' not in columns) | ('Samples' not in columns) | ('p_value' not in columns) | ('r_squared' not in columns) |('Test_Name' not in columns):
        print("ERROR:  The network file does not seem to be  KINC tidy file. It is missing one or more of the following column headers: Source, Target, Samples, p_value, r_squared or Test_Name. Please check the file.")
        exit(1)

    return net


def get_iGraph(net, directional):
    """
    Converts the KINC network dataframe into an iGraph object.

    Igraph objects are handly for performing network statistics such as
    transitivity and degree calculations.

    net :  The network dataframe created by the load_network function.

    return : An igraph object of the network loaded with the source, target and
             Similarity_Score (as the weight)
    """
    if (directional == 'Yes' or directional == 'yes'):
        g = ig.Graph(directed = True)
    else:
        g = ig.Graph()

    # Add the nodes
    v = pd.concat([net['Source'], net['Target']]).unique()
    g.add_vertices(v)

    # Add the edges
    g.add_edges(net[['Source', 'Target']].values)

    # Add any additional attributes
    for attr in net.columns:
        if attr == 'Source' or attr == 'Target':
            continue
        g.es[attr] = net[attr]

    return g

def calculate_2d_layout(g, iterations):
    """
    Calculates a typical 2D layout for the network.

    The first time this function is called on a network it may take some time
    depending on the size of the network.  The layout is saved in a file with
    the same name as the network but with a '.glayout.txt' extension in
    the working directory. On subsequent runs of this program that file is
    imported if it exists.

    g :  The network igraph object

    return : a Pandas dataframe containing the layout coordinates for
             the nodes in the network. The dataframe contains X, and Y
             dimenstional coordinates.
    """
    # Simplify the graph prior to layout calculation
    g2 = g.copy()
    g2.simplify()

    forceatlas2 = ForceAtlas2(
        # Behavior alternatives
        outboundAttractionDistribution=True,  # Dissuade hubs
        linLogMode=False,  # NOT IMPLEMENTED
        adjustSizes=False,  # Prevent overlap (NOT IMPLEMENTED)
        edgeWeightInfluence=1.0,

        # Performance
        jitterTolerance=1.0,  # Tolerance
        barnesHutOptimize=True,
        barnesHutTheta=1.2,
        multiThreaded=False,  # NOT IMPLEMENTED

        # Tuning
        scalingRatio=2.0,
        strongGravityMode=False,
        gravity=1,

        # Log
        verbose=True)

    print("Calculating 2D layout.")
    layout = pd.DataFrame(forceatlas2.forceatlas2_igraph_layout(g2, iterations=iterations).coords)
    layout.columns = ['x', 'y']
    return layout



def main():
    """
    The main function.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--net', dest='net_path', type=str, required=True, help="(required) The path to the KINC-derived network file")

    parser.add_argument('--graph_name', dest='graph_name', type=str, required=True, help="(required) The name of this graph.")

    parser.add_argument('--graph_id', dest='graph_id', type=str, required=True, help="(required) A unique name (preferrably a machine readable name) for this network.")

    parser.add_argument('--directional', dest='directional', type=str, required=False, help="Indicates if the graph is directional. Defaults to False. Set to 'Yes' if this network is directional. Default is No", default="No", choices=['Yes', 'No', 'yes', 'no']);

    parser.add_argument('--iterations', dest='iterations', type=int, default=100, help="(optional). The number of iterations to perform when calculating the Force Atlas2 layout.  This argument is only used the first time a network is viewed or if the --redo_layout argument is provided.")

    parser.add_argument('--outfile', dest='outfile', type=str, required=True, help="The name of the output graphml file.")

    args = parser.parse_args()

    # Make sure the paths exist
    if not os.path.exists(args.net_path):
        print ("ERROR: The network file cannot be found: {}".format(args.net_path))
        exit(1)

    print("Reading network file...")
    net = load_network(args.net_path)
    g = get_iGraph(net, args.directional)
    g["name"] = args.graph_name
    g["uname"] = args.graph_id

    # Add some additional vertex attributes
    g.vs['degree'] = g.degree()

    layout = calculate_2d_layout(g, args.iterations)
    g.vs['x'] = layout['x']
    g.vs['y'] = layout['y']

    g.write_graphml(args.outfile)

if __name__ == "__main__":
    main()
