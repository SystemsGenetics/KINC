#!/usr/bin/env python3
"""
Creates a Dash application that provides 3D visualization of a KINC network.

This script accepts the following arguments:

    --net : (required) The path to the KINC-derived network file

    --emx : (retuired) The path to the log2 transformed Gene Expression Matrix
            or Metabolite abundance matrix.

    --amx : (required) The path to the tab-delimited annotation matrix.
            The matrix must have at least one column that contains a unique
            list of sample names.

    --nmeta : (optional) The path to a tab-delimited node meta data file. The
              format of the file must have 4 columns, with the first
              containing the node name, the second a controlled vocabulary
              term ID, the second the term definition and the fourth the vocubulary name.

    --sample_col: (optional) The name of the column in the annotation matrix
                  that contains the unique sample names.  Defaults to "Sample"

    --debug:  (optional).  Add this argument to enable Dash application
              debugging mode.

    --iterations: (optional). The number of iterations to perform when
                  calculating the Force Atlas2 layout.  This argument is only
                  used the first time a network is viewed or if the
                  --redo_layout argument is provided.

    --redo-layout :  (optional). If the 2D and 3D network layout has already
                     been constructed it will be loaded from a file. Add this
                     arugment to force the layouts to be rebuilt and not loaded
                     from the files. To prevent Dash from rerunning the layout
                     on callbacks, this option results in the program
                     terminating. To view the application, restart without
                     this option.

"""

import argparse
import numpy as np
import pandas as pd
import igraph as ig
import plotly as py
import seaborn as sns
import plotly.graph_objects as go
from fa2 import ForceAtlas2
import random
import dash
import dash_core_components as dcc
import dash_html_components as html
import os
import json
import re
import ast
import time
from progress.bar import IncrementalBar





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





def load_gem(file_path):
    """
    Imports the tab-delimited Gene Expression Matrix (GEM) or Metabolite

    GEM files can be generated from RNA-seq data using GEMmaker. Alternatively,
    this can be a metabolite abundance matrix.

    file_path : The path to the GEM file.  The file should be log2 transformed.

    return : A pandas dataframe containing the GEM.
    """

    gem = pd.read_csv(file_path, sep="\t")
    return gem





def load_amx(file_path, sample_col = 'Sample'):
    """
    Imports the tab-delimited annotation matrix (amx).

    The matrix must have at least one column that contains a unique list of
    sample names.

    file_path :  The path to the annotation matrix.

    sample_col : The name of the column that contains the sample names. Defaults
                 to 'Sample'

    return : A pandas dataframe containing the annotation matrix.
    """
    amx = pd.read_csv(file_path, sep="\t")
    amx.index = amx[sample_col]
    return amx





def load_node_meta(file_path):
    """
    Imports the tab-delimited node metadata file.

    The format of the file must have 4 columns, with the first containing the
    node name, the second a controlled vocabulary term ID, the second the
    term definition and the fourth the vocubulary name.
    """
    nmeta = pd.read_csv(file_path, sep="\t")
    nmeta.columns = ['Node', 'Term', 'Definition', 'Vocabulary']
    nmeta.index = nmeta['Node']
    return nmeta





def get_iGraph(net):
    """
    Converts the KINC network dataframe into an iGraph object.

    Igraph objects are handly for performing network statistics such as
    transitivity and degree calculations.

    net :  The network dataframe created by the load_network function.

    return : An igraph object of the network loaded with the source, target and
             Similarity_Score (as the weight)
    """
    g = ig.Graph()

    # Add the nodes
    v = pd.concat([net['Source'], net['Target']]).unique()
    g.add_vertices(v)

    # Add the edges
    g.add_edges(net[['Source', 'Target']].values)

    # Add the edge w
    #g.es['weight'] = net['Similarity_Score']

    return g





def calculate_2d_layout(net, net_prefix, redo_layout, iterations):
    """
    Calculates a typical 2D layout for the network.

    The first time this function is called on a network it may take some time
    depending on the size of the network.  The layout is saved in a file with
    the same name as the network but with a '.glayout.txt' extension in
    the working directory. On subsequent runs of this program that file is
    imported if it exists.

    net :  The network dataframe created by the load_network function.

    net_prefix :  The filename of the file that will house the layout
                  after it is calculated. The file will be saved with this name
                  and the extension ".2Dlayout.txt"

    redo_layout :  A boolean indicting if the layout should be rebuilt rather
                   than loading from file if one exists already.

    return : a Pandas dataframe containing the layout coordinates for
             the nodes in the network. The dataframe contains X, and Y
             dimenstional coordinates.
    """

    g = get_iGraph(net)
    t = pd.Series(g.transitivity_local_undirected(), index=g.vs['name'])
    d = pd.DataFrame(g.degree(), index=g.vs['name'], columns=['Degree'])

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

    if (redo_layout | (not os.path.exists(net_prefix + '.2Dlayout.txt'))):
        print("Calculating 2D layout.")
        glayout = pd.DataFrame(forceatlas2.forceatlas2_igraph_layout(g, iterations=iterations).coords)
        glayout.columns = ['X', 'Y']
        glayout.index = g.vs['name']
        glayout = pd.concat([glayout, d, t], axis=1, sort=False)
        glayout.columns = ['X', 'Y', 'Degree', 'CC']
        glayout.to_csv(net_prefix + '.2Dlayout.txt')
    else:
        glayout = pd.read_csv(net_prefix + '.2Dlayout.txt', index_col=0)

    return glayout





def bin_edges(net):
    """
    Calculates a set of bins using the Similarity score and P-value.

    It is from these bins that the edges and nodes of the network will be
    stacked in the z-axis of the 3D plot and or colored.  Four new
    columns are added to the provided network: 'Edge_Bin', 'Pval_Bin',
    'Rsqr_Bin' and 'Relationship'.

    net :  The network dataframe created by the load_network function.

    """
    net['Edge_Bin'] = np.around(np.abs(net['Similarity_Score']), decimals=2)
    net['Pval_Bin'] = np.round(-np.log10(net['p_value']))
    net['Rsqr_Bin'] = np.around(np.abs(net['r_squared']), decimals=1)
    net['Relationship'] = np.ceil(net['Similarity_Score']).astype('str')
    net['Relationship'] = net['Relationship'].replace("-0.0", 'Negative')
    net['Relationship'] = net['Relationship'].replace("1.0", 'Positive')





def get_vertex_zlayers(net, glayout, net_prefix, redo_layout):
    """
    Uses the 2D layout and calculates the Z-coordinate for the nodes.

    net :  The network dataframe created by the load_network function.

    glayout : The dataframe containing the 2D layout of the nodes.

    net_prefix:  The filename of the file that will house the vertex layout
                 after it is calculated. The file will be saved with this name
                 and the extension ".3Dvlayers.txt"

    redo_layout :  A boolean indicting if the layout should be rebuilt rather
                   than loading from file if one exists already.

    return : A Pandas dataframe containing the X, Y and Z coordinates for the
             nodes as well as the Degree and CC (clustering coefficient) for
             each node.
    """

    def find_vlayers(row, vtype='Source', bar=None):
        if bar:
            bar.next()
        node = glayout.loc[row[vtype]]
        ebin = row['Edge_Bin']
        pbin = row['Pval_Bin']
        rbin = row['Rsqr_Bin']
        rel = row['Relationship']
        test = row['Test_Name']
        return(row[vtype], node['X'], node['Y'], ebin, pbin, rbin, rel, test, node['Degree'], node['CC'])


    if (redo_layout | (not os.path.exists(net_prefix + '.3Dvlayers.txt'))):
        print("Calculating 3D vertex layout.")
        bar = IncrementalBar('', max=net.shape[0]*2, suffix='%(percent)d%%')
        lsource = net.apply(find_vlayers, vtype='Source', bar=bar, axis=1)
        ltarget = net.apply(find_vlayers, vtype='Target', bar=bar, axis=1)
        print("")

        columns = ['Vertex', 'X', 'Y', 'EBin', 'PBin', 'RBin', 'Rel', 'Test_Name', 'Degree', 'CC']
        vlayers = pd.DataFrame.from_records(lsource.append(ltarget).values, columns=columns)
        vlayers = vlayers[vlayers.duplicated() == False]
        # We want to place the node in the layer where it first appears.
        vlayers = vlayers.groupby(by=['Vertex']).apply(lambda g: g[g['EBin'] == g['EBin'].max()])
        vlayers.reset_index(inplace=True, drop=True)
        vlayers.to_csv(net_prefix + '.3Dvlayers.txt')

    else:
        vlayers = pd.read_csv(net_prefix + '.3Dvlayers.txt', index_col=0)

    return vlayers





def get_edge_zlayers(net, glayout, net_prefix, redo_layout):
    """
    Uses the 2D layout and calculates the Z-coordinate for the edges.

    Edges are drawn as lines in the 3D scatterplot, therefore this function
    calculates the start and stop coordinates for the edges in the format
    required by the scatter3d viewer.

    net :  The network dataframe created by the load_network function.

    glayout : The dataframe containing the 2D layout of the nodes.

    net_prefix:  The filename of the file that will house the vertex layout
                 after it is calculated. The file will be saved with this name
                 and the extension ".3Delayers.txt"

    redo_layout :  A boolean indicting if the layout should be rebuilt rather
                   than loading from file if one exists already.

    return : A Pandas dataframe containing the X, Y and Z coordinates arrays
             for the edges as well as Source, Target and Samples values from
             the original network.  The X, Y and Z coordiantes are tuples.
    """

    def place_elayers(row, bar = None):
        if bar:
            bar.next()
        ebin = row['Edge_Bin']
        pbin = row['Pval_Bin']
        rbin = row['Rsqr_Bin']
        rel = row['Relationship']
        test = row['Test_Name']
        source = glayout.loc[row["Source"]]
        target = glayout.loc[row["Target"]]
        return([[source['X'], target['X'], None],
                [source['Y'], target['Y'], None],
                row["Source"],
                row["Target"],
                row["Samples"],
                ebin, pbin, rbin, rel, test])

    if (redo_layout | (not os.path.exists(net_prefix + '.3Delayers.txt'))):
        print("Calculating 3D edge layout.")
        bar = IncrementalBar('', max=net.shape[0], suffix='%(percent)d%%')
        ledge = net.apply(place_elayers, bar=bar, axis=1)
        print("")

        elayers = pd.DataFrame.from_records(ledge, columns=['X', 'Y', 'Source', 'Target', 'Samples', 'EBin', 'PBin', 'RBin', 'Rel', 'Test_Name'])
        elayers['name'] = elayers['Source'] + " (co) " + elayers['Target']
        elayers.to_csv(net_prefix + '.3Delayers.txt')
    else:
        elayers = pd.read_csv(net_prefix + '.3Delayers.txt', index_col=0)
        elayers['X'] = elayers['X'].apply(ast.literal_eval)
        elayers['Y'] = elayers['Y'].apply(ast.literal_eval)

    return elayers





def create_network_plot(net, vlayers, elayers, color_by = 'Score', layer_by = 'Score',
       camera = None, aspect = None):
    """
    Uses Plotly to create the interactive 3D visualization of the network.

    This function uses the Scatter3D plot to draw the network.  The axes are
    hidden so it appears as a typical network view.  It defaults to
    a straight on view as the network would be seen in a typical 2D viewer like
    Cytoscape.

    net :  The network dataframe created by the load_network function.

    vlayers : The dataframe containing the 3D coordinates for the nodes.

    elayers : The dataframe containing the 3D coordinates for the edges.

    camera : A dictionary containing the figure camera coordinates.

    return : a Plotly figure object.
    """

    # Default Z-indexs for lines/points to the Score value.
    Z = vlayers['EBin']
    if layer_by == 'Score':
        Z = vlayers['EBin']
    if layer_by == 'P-value':
        Z = vlayers['PBin']
    if layer_by == 'R^2':
        Z = vlayers['RBin']
    if layer_by == 'Test Name':
        Z = vlayers['Test_Name']
    if layer_by == 'Relationship':
        Z = vlayers['Rel']
    # Add the network nodes as the first trace.
    fig1 = go.Figure(data=[go.Scatter3d(x=vlayers['X'], y=vlayers['Y'],
                   z=Z, mode='markers',
                   opacity = 0.5,
                   marker=dict(symbol='circle', size=np.log10(vlayers['Degree'])*4,
                               line=dict(width=1, color="#888888")),
                   text="Node: " + vlayers['Vertex'],
                   hoverinfo='text', name='Nodes')])

    # Add the edges and bin them
    include_slider = True
    if color_by == 'Score':
        slider_title = 'Similarity Score'
    if color_by == 'P-value':
        slider_title = '-log10(p)'
    if color_by == 'R^2':
        slider_title = 'R-squared'
    if color_by == 'Test Name':
        slider_title = 'Test Name'
        include_slider = False
    if color_by == 'Relationship':
        slider_title = 'Relationship Type'
        include_slider = False

    layer_title = layer_by
    if layer_by == 'P-value':
        layer_title = '-log10(p)'

    (colorway, sliders, nticks) = create_binned_network_figure(fig1, elayers, color_by,
                            layer_by, slider_title, include_slider)

    fig1.update_layout(
        height=600,
        title=dict(text = "3D Network View", font = dict(color='#FFFFFF')),
        showlegend=True,
        legend=dict(font = dict(color="#FFFFFF")),
        margin=dict(l=10, r=10, t=30, b=10),
        paper_bgcolor="#000000",
        colorway=colorway,
        scene=dict(
          aspectmode="cube",
          xaxis=dict(showbackground=False, showline=False, zeroline=False, showgrid=False,
                     showticklabels=False, title='', showspikes=False),
          yaxis=dict(showbackground=False, showline=False, zeroline=False, showgrid=False,
                     showticklabels=False, title='', showspikes=False),
          zaxis=dict(showbackground=False, showline=False, zeroline=False, showgrid=False,
                     showticklabels=True, tickmode="auto", nticks=nticks, title=layer_title, showspikes=False, color="#FFFFFF")
        ),
        hovermode='closest',
        annotations=[dict(showarrow=False, text="", xref='paper', yref='paper',
                        x=0, y=0.1, xanchor='left', yanchor='bottom', font=dict(size=14))
                    ],
        sliders=sliders,
        autosize=True,
    )

    # We want an orthographic layout so that when looking above the edges line up
    # with the nodes.
    fig1.layout.scene.camera.projection.type = "orthographic"
    fig1.layout.scene.camera.eye = dict(x=0, y=0, z=2)
    if camera:
        fig1.layout.scene.camera.eye = camera['eye']

    fig1.layout.scene.aspectmode = 'manual'
    if aspect:
        fig1.layout.scene.aspectratio = aspect

    return fig1





def create_binned_network_figure(figure, elayers, color_by = 'Score',
    layer_by = 'Score', slider_title = 'Similarity Score', include_slider = True):

    """
    Adds the traces for the network figure based on the bin column.

    """

    color_col = 'EBin'
    if color_col == 'Score':
        color_col = 'EBin'
    if color_by == 'P-value':
        color_col = 'PBin'
    if color_by == 'R^2':
        color_col = 'RBin'
    if color_by == 'Test Name':
        color_col = 'Test_Name'
    if color_by == 'Relationship':
        color_col = 'Rel'

    layer_col = 'Edge_Bin'
    if layer_by == 'Score':
        layer_col = 'EBin'
    if layer_by == 'P-value':
        layer_col = 'PBin'
    if layer_by == 'R^2':
        layer_col = 'RBin'
    if layer_by == 'Test Name':
        layer_col = 'Test_Name'
    if layer_by == 'Relationship':
        layer_col = 'Rel'

    # Add edge traces to the figure, one each per bin.
    layer_bins = np.flip(np.sort(elayers[layer_col].unique()))
    color_bins = np.flip(np.sort(elayers[color_col].unique()))
    for bin in color_bins:
        if (not type(bin) == str):
            if (bin.dtype == "float64") & (np.isnan(bin)):
                continue

        bin_edges = elayers[elayers[color_col] == bin]

        # Reformat the elayers for use by the Scatter3d function.
        eX = np.hstack(bin_edges['X'])
        eY = np.hstack(bin_edges['Y'])
        eZ = np.hstack(bin_edges[layer_col].repeat(3))
        names = bin_edges['name'][bin_edges.index.repeat(3)]

        # Create the scatterplot containing the lines for edges.
        figure.add_trace(go.Scatter3d(x=eX, y=eY, z=eZ,
                       mode='lines',
                       line=dict(width=1),
                       text="Edge: " + names,
                       hoverinfo='text', name=bin,
                       customdata=bin_edges.index.repeat(3)))

    #  Add a slider for the network viewer
    if include_slider:
        steps = []
        steps.append(dict(
            method="restyle",
            args=["visible", [True] * (len(color_bins) + 2)],
            label='all'
        ))
        steps.append(dict(
            method="restyle",
            args=["visible", [False] * (len(color_bins) + 2)],
            label='nodes'
        ))
        steps[1]["args"][1][0] = True
        for i in range(len(color_bins)):
            step = dict(
                method="restyle",
                args=["visible", [False] * (len(color_bins) + 2)],
                label=color_bins[i]
            )
            # Turn on the layers for this step and leave on the nodes layer.
            step["args"][1][0] = True
            for j in range(1,i+2):
                step["args"][1][j] = True

            # Set the label.
            steps.append(step)


        colorway = ["#FFFFFF"] + sns.color_palette('viridis_r', color_bins.size).as_hex()

        sliders = [dict(
            active=0,
            currentvalue={"prefix": slider_title + ": "},
            pad={"t": 50},
            steps=steps,
            font=dict(color = '#FFFFFF'),
            tickcolor='#FFFFFF',
            len=1)]
    else:
        colorway = ["#FFFFFF"] + sns.color_palette('muted', color_bins.size).as_hex()
        sliders = []

    nticks = layer_bins.size
    if layer_by == 'Score':
        nticks = int(nticks / 2)
    return (colorway, sliders, nticks)




def create_expression_scatterplot(gem, amx, elayers, color_col=None, edge_index = None):
    """
    Uses Plotly to create the interactive 3D scatterplot of co-expression

    This function uses the Scatter3D plot to draw the co-expression scatterplot.
    It defaults to a straight on view but can be interactively rotated,
    panned, etc.

    net :  The network dataframe created by the load_network function.

    amx : The annotation matrix dataframe created by the load_amx function.

    elayers : The dataframe containing the 3D coordinates for the edges.

    color_col : The name of the column in the amx that contains the category
                that should be used for coloring the points in the plot.

    edge_index : The numerical index of the edge in the elayers dataframe
                 that is to be plotted.

    return : a Plotly figure object.
    """
    if edge_index is None:
        return go.Figure(go.Scatter3d())

    node1 = elayers.iloc[edge_index]['Source']
    node2 = elayers.iloc[edge_index]['Target']
    samples = elayers.iloc[edge_index]['Samples']


    # Generate the dataframe for the expression scatterplot
    sdata = pd.DataFrame(dict(X=gem.loc[node1].values, Y=gem.loc[node2].values))
    sdata.index = gem.columns
    sdata = sdata.join(amx, how='left')

    # Calculate the sizes of the points.
    sizes = pd.Series(list(samples))
    sizes = sizes.replace(to_replace=r'[^1]', value='8', regex=True)
    sizes = sizes.replace({'1': '16'})
    sizes = sizes.astype('int')
    sizes.index = sdata.index

    # Generate the colors for the samples.
    if (color_col == None):
        color_col = 'Cluster'

    # If the column is 'Cluster' we need to add it to the dataframe. The
    # Cluster column simply lists if the sample is in the cluster or not.
    if (color_col == 'Cluster'):
        inout = pd.Series(list(samples))
        inout = inout.replace(to_replace=r'[^1]', value='Out', regex=True)
        inout = inout.replace({'1': 'In'})
        inout.index = gem.columns
        sdata = pd.concat([sdata, inout.rename('Cluster')], 1)

    # Is this a categorical column?
    is_categorical = False
    categories = sdata[color_col].unique()
    if (categories.dtype == object):
        is_categorical = True

    # Now draw the plot
    nticks = None
    tickmode = 'auto'
    ticktext = None
    tickvals = None
    if is_categorical:
        num_categories = categories.shape[0]
        tickmode = 'array'
        ticktext = categories
        tickvals = np.arange(0, num_categories) / (num_categories - 1) - 0.5
        replace_df = pd.DataFrame({'Categories' : categories,'Z' : tickvals})
        sdata['Z'] = sdata[color_col].replace(
                             to_replace=replace_df['Categories'].values,
                             value=replace_df['Z'].values)

        nticks = num_categories
        showlegend = True
        first_category = (sdata[color_col] == categories[0])
        fig2 = go.Figure(data=[go.Scatter3d(x=sdata[first_category]['X'],
                        z=sdata[first_category]['Y'],y=sdata[first_category]['Z'],
                        mode='markers',
                        marker=dict(symbol='circle', size=sizes[first_category]),
                        text= sdata[first_category].index, hoverinfo='text',
                        name=str(categories[0]))])

        for i in range(1, len(categories)):
            next_category = (sdata[color_col] == categories[i])
            fig2.add_trace(go.Scatter3d(x=sdata[next_category]['X'],
                            z=sdata[next_category]['Y'], y=sdata[next_category]['Z'],
                            mode='markers',
                            marker=dict(symbol='circle',size=sizes[next_category]),
                            text= sdata[next_category].index,
                            hoverinfo='text', name=str(categories[i])))
    else:
        num_categories = None
        sdata['Z'] = sdata[color_col]
        tickvals = []
        showlegend = False
        fig2 = go.Figure(data=[go.Scatter3d(x=sdata['X'], z=sdata['Y'], y=sdata['Z'],
                        mode='markers',
                        marker=dict(symbol='circle', size=sizes,
                                    color=sdata['Z'], colorscale='Viridis'),
                        text= sdata.index, hoverinfo='text')])

    fig2.update_layout(
        height=600,
        title="3D Edge Co-Expression Scatterplot",
        showlegend=showlegend,
        legend={'itemsizing': 'constant'},
        margin=dict(l=10, r=10, t=30, b=10),
        scene=dict(
          aspectmode="cube",
          xaxis=dict(showbackground=True, showline=True, zeroline=True, showgrid=True,
                     showticklabels=True, title=node1,
                     showspikes=True),
          zaxis=dict(showbackground=True, showline=True, zeroline=True, showgrid=True,
                     showticklabels=True, title=node2,
                     showspikes=True),
          yaxis=dict(showbackground=True, showline=True, zeroline=True, showgrid=True,
                     showticklabels=True, title=color_col,
                     tickmode=tickmode, ticktext=ticktext, tickvals=tickvals, nticks=nticks, showspikes=True),
        ),
        hovermode='closest',
        annotations=[dict(showarrow=False,
                        text="",
                        xref='paper', yref='paper',
                        x=0, y=0.1, xanchor='left', yanchor='bottom', font=dict(size=14))
                    ],
        datarevision = time.time()
    )

    fig2.layout.scene.camera.projection.type = "orthographic"
    fig2.layout.scene.camera.eye = dict(x=0, y=-1, z=0)

    return fig2




def create_dash_edge_table(net, edge_index = None):
    """
    Constructs the HTML table that holds edge information for the Dash appself.

    elayers : The dataframe containing the 3D coordinates for the edges.

    edge_index : The numerical index of the edge in the elayers dataframe
                 that is to be plotted.

    returns : a Dash html.Table object.
    """

    htr_style = {'background-color' : '#f2f2f2'}
    hth_style = {'text-align' : 'left',
                 'padding' : '5px',
                 'border-bottom' : '1px solid #ddd'};
    th_style = {'text-align' : 'left',
                'padding' : '5px',
                'border-bottom' : '1px solid #ddd'};

    net_fixed = net.drop(['Samples', 'Edge_Bin', 'Pval_Bin', 'Rsqr_Bin', 'Relationship'], axis=1)
    if ('p_value' in net_fixed.columns):
        net_fixed['p_value'] = net_fixed['p_value'].apply(np.format_float_scientific, precision=4)

    columns = net_fixed.columns
    table_rows = []
    table_rows.append(html.Tr([html.Th(col, style=hth_style) for col in columns], style=htr_style))
    if not edge_index == None:
        row_vals = net_fixed.iloc[edge_index]
        source = row_vals['Source']
        target = row_vals['Target']
        row_vals = net_fixed[(net_fixed['Source'] == source) & (net_fixed['Target'] == target)]
        for index, row in row_vals.iterrows():
            table_rows.append(html.Tr([html.Th(row[col], style=th_style) for col in columns]))

    return html.Table(children=table_rows)





def create_dash_sample_table(net, amx, sample = None):
    """
    Constructs the HTML table that holds sample information for the Dash appself.

    amx : The annotation matrix dataframe created by the load_amx function.

    sample : The name of the sample to display

    returns : a Dash html.Table object.
    """

    htr_style = {'background-color' : '#f2f2f2'}
    hth_style = {'text-align' : 'left',
                 'padding' : '5px',
                 'border-bottom' : '1px solid #ddd'};
    th_style = {'text-align' : 'left',
                'padding' : '5px',
                'border-bottom' : '1px solid #ddd'};

    columns = amx.columns
    table_rows = []
    table_rows.append(html.Tr([html.Th(col, style=hth_style) for col in columns], style=htr_style))
    if sample:
        row_vals = amx.loc[sample]
        table_rows.append(html.Tr([html.Th(row_vals[col], style=th_style) for col in columns]))

    return html.Table(children=table_rows)





def create_dash_node_table(net, nmeta, nodes = None):
    """
    Constructs the HTML table that holds node information for the Dash app.

    net :  The network dataframe created by the load_network function.

    nmeta : The dataframe containing the node metadata.

    node : The name of the node to display

    returns : a Dash html.Table object.
    """

    htr_style = {'background-color' : '#f2f2f2'}
    hth_style = {'text-align' : 'left',
                 'padding' : '5px',
                 'border-bottom' : '1px solid #ddd'};
    th_style = {'text-align' : 'left',
                'padding' : '5px',
                'border-bottom' : '1px solid #ddd'};


    if nmeta:
        columns = nmeta.columns
        table_rows = []
        table_rows.append(html.Tr([html.Th(col, style=hth_style) for col in columns], style=htr_style))
        if (not nmeta is None) & (not nodes is None):
            row_vals = nmeta.loc[nodes]
            for index, row in row_vals.iterrows():
                table_rows.append(html.Tr([html.Th(row[col], style=th_style) for col in columns]))

        return html.Table(children=table_rows)
    return ''






def create_condition_select(amx, sample_col = 'Cluster'):
    """
    Creates a Dash select dropdown for selecting the condition to view.

    This dropbox is intended to change the 3D co-expression scatterplot.

    amx : The annotation matrix dataframe created by the load_amx function.

    color_col : The name of the column in the amx that contains the category
                that should be used for coloring the points in the plot.

    return : A Dash dcc.Dropdown object.

    """
    columns = np.sort(amx.columns.values)

    # Holds the list of columns to keep.
    keep = []
    keep.append('Cluster')

    # Exclude any columns with just a single value or any columns with as
    # many unique values as there are elements
    for col in columns:
        if len(amx[col].dropna().unique()) <= 1:
            continue
        if len(amx[col].dropna().unique()) == amx[col].size:
            continue
        keep.append(col)

    # Build the select element.
    select = dcc.Dropdown(
        id = 'coexp-condition-select',
        options = [
          {'label' : col, 'value' : col} for col in keep
        ],
        value = 'Cluster')
    return select





def create_edge_color_select(net):
    """
    Creates a Dash select dropdown for selecting the network attribute to view.

    This dropbox is intended to change the 3D network layout view.

    net :  The network dataframe created by the load_network function.

    return : A Dash dcc.Dropdown object.

    """

    options = ['Score']
    if 'p_value' in net.columns:
        options.append('P-value')
    if 'Test_Name' in net.columns:
        options.append('Test Name')
    if 'r_squared' in net.columns:
        options.append('R^2')
    options.append('Relationship')

    select = dcc.Dropdown(
        id = 'edge-color-select',
        options = [
            {'label' : col, 'value' : col} for col in options
        ],
        value = 'Score'
    )
    return select






def create_edge_layer_select(net):
    """
    Creates a Dash select dropdown for selecting the network attribute to view.

    This dropbox is intended to change the 3D network layout view.

    net :  The network dataframe created by the load_network function.

    return : A Dash dcc.Dropdown object.

    """

    options = ['Score']
    if 'p_value' in net.columns:
        options.append('P-value')
    if 'Test_Name' in net.columns:
        options.append('Test Name')
    if 'r_squared' in net.columns:
        options.append('R^2')
    options.append('Relationship')

    select = dcc.Dropdown(
        id = 'edge-layer-select',
        options = [
            {'label' : col, 'value' : col} for col in options
        ],
        value = 'Score'
    )
    return select



def build_application(net, gem, amx, nmeta, vlayers, elayers, sample_col,
    net_name):

    """
    Creates the Dash application.

    The Dash application will provide all of the interactive plots, tables and
    filters to interacitvely exploring the network.

    net :  The network dataframe created by the load_network function.

    gem :  The GEM dataframe created by the load_gem function.

    amx : The annotation matrix dataframe created by the load_amx function.

    nmeta : The dataframe containing the node metadata.

    vlayers : The dataframe containing the 3D coordinates for the nodes.

    elayers : The dataframe containing the 3D coordinates for the edges.

    sample_col : The name of the column in the amx that contains the sample
                 name.

    net_name : The name of the network to display.


    return : The Dash application object.
    """
    app = dash.Dash()
    app.layout = html.Div([
        # Header Row
        html.Div(className='row', id = "header", children=[
            html.Img(
                src="https://raw.githubusercontent.com/SystemsGenetics/KINC/master/docs/images/kinc.png",
                style={"height" : "55px","display" : "inline-block",
                  "padding" : "0px", "margin" : "0px 10px 0px 10px"}),
            html.H1(children="3D Network Explorer",
                style={
                  "display" : "inline-block", "padding" : "10px 0px 0px 0px",
                  "margin" : "0px", "vertical-align" : "top"}),
            html.Div(children="Network name: " + net_name,
                style={"padding" : "0px 0px 0px 10px"}),
        ]),
        # Graph Row
        html.Div(children=[
            # 3D network plot
            html.Div(children=[
                html.Div(id='edge-select-box', children=[
                    html.Label('Color', style={'color': '#FFFFFF'}),
                    create_edge_color_select(net),
                    html.Label('Layer', style={'color': '#FFFFFF'}),
                    create_edge_layer_select(net)],
                    style={'padding-bottom' : '10px'}),
                dcc.Graph(id = 'network-3dview',
                  figure = create_network_plot(net, vlayers, elayers),
                  config =  {
                        'toImageButtonOptions' : {
                          'filename': 'kinc_3d_network_view',
                          'width': 800,
                          'height': 600,
                          'format': 'svg',
                          'scale' : 2
                         }
                       }),
                dcc.Input(
                    id='current-network-3dview-camera', type="number",
                    value=0, style= {'display' : 'none'}),
                dcc.Input(
                    id='current-network-3dview-aspect', type="number",
                    value=0, style= {'display' : 'none'})],
                style = {"width" : "55%", "display" : "inline-block",
                         "border" : "1px solid black", "padding" : "10px",
                         "background-color" : "black", "margin" : "10px"},
            ),
            # 3D Co-Expression scatterplot.
            html.Div(children=[
                html.Div(id='coexp-condition-select-box', children=[
                    create_condition_select(amx, sample_col)],
                    style={'padding-bottom' : '10px'}),
                dcc.Graph(id = 'edge-expression-3dview',
                    figure = create_expression_scatterplot(gem, amx, elayers),
                    config =  {
                          'toImageButtonOptions' : {
                            'filename': 'kinc_3d_expression_scatterplot',
                            'width': 800,
                            'height': 600,
                            'format': 'svg',
                            'scale' : 1
                           }
                         }),
                dcc.Input(
                    id='selected-edge', type="number",
                    value=-1, style= {'display' : 'none'}),
                dcc.Input(
                    id='current-expr-camera-coords', type="number",
                    value=0, style= {'display' : 'none'})],
                style = {"width" : "35%", "display" : "inline-block",
                         "border" : "1px solid black", "padding" : "10px",
                         "margin" : "10px", "vertical-align" : "top"},
            ),
        ]),
        # Table Row
        html.Div(className='row', children=[
            html.H3(children="Edge Details", style={
                  "display" : "inline-block", "padding" : "10px 0px 0px 0px",
                  "margin" : "0px", "vertical-align" : "top"}),
            html.Div(id = "edge-table", children=[
                create_dash_edge_table(net)]),
            html.H3(children="Node Detail", style={
                  "display" : "inline-block", "padding" : "10px 0px 0px 0px",
                  "margin" : "0px", "vertical-align" : "top"}),
            html.Div(id = "node-table", children=[
                create_dash_node_table(net, nmeta)]),
            html.H3(children="Sample Detail", style={
                  "display" : "inline-block", "padding" : "10px 0px 0px 0px",
                  "margin" : "0px", "vertical-align" : "top"}),
            html.Div(id = "sample-table", children=[
                create_dash_sample_table(net, amx)]),
            # Holds the values of the clicked node or cliced edge nodes.
            dcc.Input(
                id='selected-nodes', type="text",
                value=0, style= {'display' : 'none'}),
        ]),
    ])


    # Callback when an object in the network plot is clicked.
    @app.callback(
        dash.dependencies.Output('selected-edge', 'value'),
        [dash.dependencies.Input('network-3dview', 'clickData')])
    def set_current_edge(clickData):
        edge_index = None
        edge_nodes = None
        nodes = None
        if (clickData):
            points = clickData['points']
            efound = re.match('^Edge: (.*?) \(co\) (.*?)$', points[0]['text'])
            nfound = re.match('^Node: (.*?)$', points[0]['text'])
            if (efound):
                edge_index = points[0]['customdata']
                row_vals = elayers.iloc[edge_index]
                source = row_vals['Source']
                target = row_vals['Target']
                edge_nodes = [source, target]
                return edge_index #, json.dumps(edge_nodes))
            # elif (nfound):
            #     return (-1, json.dumps([nfound.group(1)]))
        raise dash.exceptions.PreventUpdate


    @app.callback(
         dash.dependencies.Output('node-table', 'children'),
         [dash.dependencies.Input('selected-nodes', 'value')])
    def update_node_table(nodes):
        if nodes == 0:
             raise dash.exceptions.PreventUpdate
        return create_dash_node_table(net, nmeta, json.loads(nodes))

    @app.callback(
        dash.dependencies.Output('sample-table', 'children'),
        [dash.dependencies.Input('edge-expression-3dview', 'clickData')])
    def update_sample_table(clickData):
        if (clickData):
            sample = clickData['points'][0]['text']
            return create_dash_sample_table(net, amx, sample)
        raise dash.exceptions.PreventUpdate

    @app.callback(
         dash.dependencies.Output('edge-table', 'children'),
         [dash.dependencies.Input('selected-edge', 'value')])
    def update_edge_table(current_edge):
        if (current_edge == -1):
          current_edge = None
        return create_dash_edge_table(net, current_edge)

    @app.callback(
        dash.dependencies.Output('current-network-3dview-camera', 'value'),
        [dash.dependencies.Input('network-3dview', 'relayoutData')])
    def set_current_camera(relayoutData):
        if (relayoutData):
            if 'scene.camera' in relayoutData.keys():
                camera = json.dumps(relayoutData["scene.camera"])
                return camera
        raise dash.exceptions.PreventUpdate

    @app.callback(
        dash.dependencies.Output('current-network-3dview-aspect', 'value'),
        [dash.dependencies.Input('network-3dview', 'relayoutData')])
    def set_network_aspect(relayoutData):
        if (relayoutData):
            if 'scene.aspectratio' in relayoutData.keys():
                aspect = json.dumps(relayoutData["scene.aspectratio"])
                return aspect
        raise dash.exceptions.PreventUpdate

    # Callback to update the co-expression plot when the edge changes.
    @app.callback(
        dash.dependencies.Output('edge-expression-3dview', 'figure'),
        [dash.dependencies.Input('selected-edge', 'value'),
         dash.dependencies.Input('coexp-condition-select', 'value')])
    def update_expression_plot(current_edge, color_col):
        if (current_edge == -1):
          current_edge = None
        return create_expression_scatterplot(gem, amx, elayers, color_col, current_edge)


    @app.callback(
        dash.dependencies.Output('network-3dview', 'figure'),
        [dash.dependencies.Input('edge-color-select', 'value'),
         dash.dependencies.Input('edge-layer-select', 'value')],
        [dash.dependencies.State('current-network-3dview-camera', 'value'),
         dash.dependencies.State('current-network-3dview-aspect', 'value')])
    def update_network_plot(color_by,  layer_by, camera_vals, aspect_vals):
        camera = None
        aspect = None
        if (type(camera_vals) == str):
            camera = json.loads(camera_vals)
        if (type(aspect_vals) == str):
            aspect = json.loads(aspect_vals)

        if not camera and not aspect:
            raise dash.exceptions.PreventUpdate

        return create_network_plot(net, vlayers, elayers, color_by, layer_by, camera, aspect)



    return app



def main():
    """
    The main function.

    """

    parser = argparse.ArgumentParser()
    parser.add_argument('--net', dest='net_path', type=str, required=True)
    parser.add_argument('--emx', dest='gem_path', type=str, required=True)
    parser.add_argument('--amx', dest='amx_path', type=str, required=True)
    parser.add_argument('--sample_col', dest='sample_col', type=str, required=False, default='Sample')
    parser.add_argument('--nmeta', dest='nmeta', type=str, required=False)
    parser.add_argument('--debug', dest='debug', action='store_true', default=False)
    parser.add_argument('--redo-layout', dest='redo_layout', action='store_true', default=False)
    parser.add_argument('--iterations', dest='iterations', type=int, default=100)
    args = parser.parse_args()

    # Make sure the paths exist
    if not os.path.exists(args.net_path):
        print ("ERROR: The network file cannot be found: {}".format(args.net_path))
        exit(1)
    if not os.path.exists(args.gem_path):
        print ("ERROR: The expression matrix file cannot be found: {}".format(args.gem_path))
        exit(1)
    if not os.path.exists(args.amx_path):
        print ("ERROR: The annotation matrix file cannot be found: {}".format(args.amx_path))
        exit(1)
    if not args.nmeta is None:
        if not os.path.exists(args.nmeta):
            print ("ERROR: The node metadata file cannot be found: {}".format(args.nmeta))
            exit(1)

    # Load the input data.
    print("Reading network file...")
    net = load_network(args.net_path)
    print("Reading GEM file...")
    gem = load_gem(args.gem_path)
    print("Reading experioment annotation file...")
    amx = load_amx(args.amx_path, args.sample_col)

    nmeta = None
    if not args.nmeta is None:
        print("Reading the node metadata file...")
        nmeta = load_node_meta(args.nmeta)


    # Get the filename of the network file minus the extension.
    (net_prefix, net_ext) = os.path.splitext(os.path.basename(args.net_path))

    # Calculate a 2D layout for the network
    glayout = calculate_2d_layout(net, net_prefix, args.redo_layout, args.iterations)

    # Calculate the Z-coorinate positions for the verticies and edges.
    bin_edges(net)
    vlayers = get_vertex_zlayers(net, glayout, net_prefix, args.redo_layout)
    elayers = get_edge_zlayers(net, glayout, net_prefix, args.redo_layout)


    # If the user requested we rebuild the layout then terminate so Dash
    # doesn't try to rebuild the layout on each callback.
    if args.redo_layout:
        print ("Layouts have been built. Please relaunch without the --redo-layouts option to view the app.")
        exit(0)

    # Launch the dash application
    print("Launching application...")
    app = build_application(net, gem, amx, nmeta, vlayers, elayers, args.sample_col, net_prefix)
    app.run_server(debug=args.debug)

    exit(0)



if __name__ == "__main__":
    main()
