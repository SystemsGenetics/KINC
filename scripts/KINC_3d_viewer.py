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

    --sample_col: (optional) The name of the column in the annotation matrix
                  that contains the unique sample names.  Defaults to "Sample"

    --debug:  (optional).  Add this argument to enable Dash application
              debugging mode.

"""

import argparse
import numpy as np
import pandas as pd
import igraph as ig
import plotly as py
import seaborn as sns
import plotly.graph_objects as go
import networkx as nx
import dash
import dash_core_components as dcc
import dash_html_components as html
import os
import json
import re





def load_network(file_path):
    """
    Imports the KINC-generated network file (either full or Tidy versions).

    file_path : The path to the network file.

    return : A pandas dataframe containing the network.
    """

    net = pd.read_csv(file_path, sep="\t")
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
    g.es['weight'] = net['Similarity_Score']

    return g





def calculate_2d_layout(net, net_prefix):
    """
    Calculates a typical 2D layout for the network.

    The first time this function is called on a network it may take some time
    depending on the size of the network.  The layout is saved in a file with
    the same name as the network but with a '.glayout.txt' extension in
    the working directory. On subsequent runs of this program that file is
    imported if it exists.

    net :  The network dataframe created by the load_network function.

    net_prefix:  The filename of the file that will house the layout
                 after it is calculated. The file will be saved with this name
                 and the extension ".2Dlayout.txt"

    return : a Pandas dataframe containing the layout coordinates for
             the nodes in the network. The dataframe contains X, and Y
             dimenstional coordinates.
    """

    g = get_iGraph(net)
    t = pd.Series(g.transitivity_local_undirected(), index=g.vs['name'])
    d = pd.DataFrame(g.degree(), index=g.vs['name'], columns=['Degree'])

    G = nx.Graph()
    net['Weight'] = np.abs(net['Similarity_Score'])
    G.add_weighted_edges_from(net[['Source','Target','Weight']].values)
    if (not os.path.exists(net_prefix + '.2Dlayout.txt')):
        glayout = pd.DataFrame(nx.drawing.layout.kamada_kawai_layout(G)).transpose()
        glayout.columns = ['X', 'Y']
        glayout = pd.concat([glayout, d, t], axis=1, sort=False)
        glayout.columns = ['X', 'Y', 'Degree', 'CC']
        glayout.to_csv(net_prefix + '.2Dlayout.txt')
    else:
        glayout = pd.read_csv(net_prefix + '.2Dlayout.txt', index_col=0)

    return glayout





def bin_edges(net):
    """
    Calculates a set of bins using the Similarity score.

    It is from these bins that the edges and nodes of the network will be
    stacked in the z-axis of the 3D plot.

    net :  The network dataframe created by the load_network function.

    return : a numpy array containing the list of bin values.
    """
    return np.around(np.abs(net['Similarity_Score']), decimals=2)





def get_vertex_zlayers(net, glayout):
    """
    Uses the 2D layout and calculates the Z-coordinate for the nodes.

    net :  The network dataframe created by the load_network function.

    glayout : The dataframe containing the 2D layout of the nodes.

    return : A Pandas dataframe containing the X, Y and Z coordinates for the
             nodes as well as the Degree and CC (clustering coefficient) for
             each node.
    """
    net['Bin'] = bin_edges(net)
    net['Weight'] = np.abs(net['Similarity_Score'])

    def find_vlayers(row, vtype='Source'):
        node = glayout.loc[row[vtype]]
        sbin = row['Bin']
        return(row[vtype], node['X'], node['Y'], sbin, node['Degree'], node['CC'])

    lsource = net.apply(find_vlayers, vtype='Source', axis=1)
    ltarget = net.apply(find_vlayers, vtype='Target', axis=1)

    vlayers = pd.DataFrame.from_records(lsource.append(ltarget).values, columns=['Vertex', 'X', 'Y', 'Z', 'Degree', 'CC'])
    vlayers.dropna(inplace=True)
    vlayers = vlayers[vlayers.duplicated() == False]
    vlayers.head()

    vlayers = vlayers.groupby(by=['Vertex']).apply(lambda g: g[g['Z'] == g['Z'].max()])
    return vlayers





def get_edge_zlayers(net, glayout):
    """
    Uses the 2D layout and calculates the Z-coordinate for the edges.

    Edges are drawn as lines in the 3D scatterplot, therefore this function
    calculates the start and stop coordinates for the edges in the format
    required by the scatter3d viewer.

    net :  The network dataframe created by the load_network function.

    glayout : The dataframe containing the 2D layout of the nodes.

    return : A Pandas dataframe containing the X, Y and Z coordinates arrays
             for the edges as well as Source, Target and Samples values from
             the original network.  The X, Y and Z coordiantes are tuples.
    """

    def place_elayers(row):
        sbin = row['Bin']
        source = glayout.loc[row["Source"]]
        target = glayout.loc[row["Target"]]
        return([[source['X'], target['X'], None],
                [source['Y'], target['Y'], None],
                [sbin, sbin, None],
                row["Source"],
                row["Target"],
                row["Samples"],
                sbin])

    ledge = net.apply(place_elayers, axis=1)

    elayers = pd.DataFrame.from_records(ledge, columns=['X', 'Y', 'Z', 'Source', 'Target', 'Samples', 'Bin'])
    elayers.dropna(inplace=True)
    elayers['name'] = elayers['Source'] + " (co) " + elayers['Target']

    return elayers





def create_network_plot(net, vlayers, elayers):
    """
    Uses Plotly to create the interactive 3D visualization of the network.

    This function uses the Scatter3D plot to draw the network.  The axes are
    hidden so it appears as a typical network view.  It defaults to
    a straight on view as the network would be seen in a typical 2D viewer like
    Cytoscape.

    net :  The network dataframe created by the load_network function.

    vlayers : The dataframe containing the 3D coordinates for the nodes.

    elayers : The dataframe containing the 3D coordinates for the edges.

    return : a Plotly figure object.
    """

    fig1 = go.Figure(data=[go.Scatter3d(x=vlayers['X'], y=vlayers['Y'],
                   z=vlayers['Z'], mode='markers',
                   marker=dict(symbol='circle', size=np.log10(vlayers['Degree'])*4),
                   text="Node: " + vlayers['Vertex'],
                   hoverinfo='text', name='Nodes')])

    # Calculate a color per layer
    bins = np.flip(np.sort(elayers['Bin'].unique()))

    # Add edge traces to the figure, one each per bin.
    for bin in bins:
        bin_edges = elayers[elayers['Bin'] == bin]

        # Reformat the elayers for use by the Scatter3d function.
        eX = np.hstack(bin_edges['X'])
        eY = np.hstack(bin_edges['Y'])
        eZ = np.hstack(bin_edges['Z'])
        names = bin_edges['name'][bin_edges.index.repeat(3)]

        # Create the scatterplot containing the lines for edges.
        fig1.add_trace(go.Scatter3d(x=eX, y=eY, z=eZ,
                       mode='lines', line=dict(width=1),
                       text="Edge: " + names,
                       hoverinfo='text', name=bin,
                       customdata=bin_edges.index.repeat(3)))

    #  Add a slider for the network viewer
    steps = []
    num_layers = len(fig1.data)
    for i in range(num_layers):
        step = dict(
            method="restyle",
            # Disable all layers.
            args=["visible", [False] * num_layers if (i > 0) else [True] * num_layers],
            label='all' if (i == 0) else bins[i-1]
        )
        # Turn on the layers for this step. It should always be the first
        # layer and i'th layer.
        step["args"][1][0] = True
        for j in range(i):
            step["args"][1][j] = True

        # Set the label.
        steps.append(step)

    fig1.update_layout(
        height=600,
        title=dict(text = "3D Network View", font = dict(color='#FFFFFF')),
        showlegend=True,
        legend=dict(font = dict(color="#FFFFFF")),
        margin=dict(l=10, r=10, t=30, b=10),
        paper_bgcolor="#000000",
        colorway=["#FFFFFF"] + sns.color_palette('viridis_r', bins.size).as_hex(),
        scene=dict(
          aspectmode="cube",
          xaxis=dict(showbackground=False, showline=False, zeroline=False, showgrid=False,
                     showticklabels=False, title='', showspikes=False),
          yaxis=dict(showbackground=False, showline=False, zeroline=False, showgrid=False,
                     showticklabels=False, title='', showspikes=False),
          zaxis=dict(showbackground=False, showline=False, zeroline=False, showgrid=False,
                     showticklabels=False, title='', showspikes=False, color="#FFFFFF")
        ),
        hovermode='closest',
        annotations=[dict(showarrow=False,
                        text="",
                        xref='paper', yref='paper',
                        x=0, y=0.1, xanchor='left', yanchor='bottom', font=dict(size=14))
                    ],
        sliders=[dict(
            active=0,
            currentvalue={"prefix": "Similarity: "},
            pad={"t": 50},
            steps=steps,
            font=dict(color = '#FFFFFF'),
            tickcolor='#FFFFFF',
            len=1)],
    )

    # We want an orthographic layout so that when looking above the edges line up
    # with the nodes.
    fig1.layout.scene.camera.projection.type = "orthographic"
    fig1.layout.scene.camera.eye = dict(x=0., y=0., z=2)

    return fig1





def create_expression_scatterplot(gem, amx, elayers, color_col=None, edge_index=0):
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
    node1 = elayers.iloc[edge_index]['Source']
    node2 = elayers.iloc[edge_index]['Target']
    samples = elayers.iloc[edge_index]['Samples']

    # Generate the dataframe for the expression scatterplot
    sdata = pd.DataFrame(dict(x=gem.loc[node1].values,y=gem.loc[node2].values))
    sdata.index = gem.columns
    sdata = sdata.join(amx, how='left')

    if (color_col == None):
        color_col = amx.columns[1]

    # Generate the colors for the samples.
    categories = sdata[color_col].unique()
    num_categories = categories.shape[0]
    color_rep = pd.DataFrame({
        'Categories' : categories,
        'Color' : np.arange(0, num_categories)})
    sdata['color'] = sdata[color_col].replace(
                         to_replace=color_rep['Categories'].values,
                         value=color_rep['Color'].values)

    # Generate the Z-coordinates for the samples.
    ticks = np.arange(0, num_categories) / (num_categories - 1) - 0.5
    color_rep = pd.DataFrame({
        'Categories' : categories,
        'Z' : ticks})
    sdata['z'] = sdata[color_col].replace(
                         to_replace=color_rep['Categories'].values,
                         value=color_rep['Z'].values)

    # Calculate the sizes of the points.
    sizes = pd.Series(list(samples))
    sizes = sizes.replace(to_replace=r'[^1]', value='8', regex=True)
    sizes = sizes.replace({'1': '16'})
    sizes = sizes.astype('int')
    sizes.index = sdata.index

    # Add the first category as trace #1.
    first_category = (sdata[color_col] == categories[0])
    fig2 = go.Figure(data=[go.Scatter3d(
                    x=sdata[first_category]['x'],
                    z=sdata[first_category]['y'],
                    y=sdata[first_category]['z'],
                    mode='markers',
                    marker=dict(symbol='circle',size=sizes[first_category]),
                    text= sdata[first_category]['Sample'],
                    hoverinfo='text',
                    name=str(categories[0]))])

    for i in range(1, len(categories)):
        next_category = (sdata[color_col] == categories[i])
        fig2.add_trace(go.Scatter3d(
                        x=sdata[next_category]['x'],
                        z=sdata[next_category]['y'],
                        y=sdata[next_category]['z'],
                        mode='markers',
                        marker=dict(symbol='circle',size=sizes[next_category]),
                        text= sdata[next_category]['Sample'],
                        hoverinfo='text',
                        name=str(categories[i])))

    fig2.update_layout(
        height=600,
        title="3D Edge Co-Expression Scatterplot",
        showlegend=True,
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
                     showticklabels=True, title='Condition',
                     tickvals=ticks,
                     ticktext=categories, showspikes=True),
        ),
        hovermode='closest',
        annotations=[dict(showarrow=False,
                        text="",
                        xref='paper', yref='paper',
                        x=0, y=0.1, xanchor='left', yanchor='bottom', font=dict(size=14))
                    ],
    )

    fig2.layout.scene.camera.projection.type = "orthographic"
    fig2.layout.scene.camera.eye = dict(x=0, y=1, z=0)

    return fig2




def create_dash_edge_table(net, edge_index = 0):
    """
    Constructs the HTML table that holds edge information for the Dash appself.

    elayers : The dataframe containing the 3D coordinates for the edges.

    edge_index : The numerical index of the edge in the elayers dataframe
                 that is to be plotted.

    returns : a Dash html.Table object.
    """

    net_fixed = net.drop('Samples', axis=1)
    if ('p_value' in net_fixed.columns):
        net_fixed['p_value'] = net_fixed['p_value'].apply(np.format_float_scientific, precision=4)
    columns = net_fixed.columns
    row = net_fixed.iloc[edge_index]

    htr_style = {'background-color' : '#f2f2f2'}
    hth_style = {'text-align' : 'left',
                 'padding' : '5px',
                 'border-bottom' : '1px solid #ddd'};
    th_style = {'text-align' : 'left',
                'padding' : '5px',
                'border-bottom' : '1px solid #ddd'};
    table = html.Table(
        children=[
            html.Tr([html.Th(col, style=hth_style) for col in columns], style=htr_style),
            html.Tr([html.Th(row[col], style=th_style) for col in columns]),
        ],
    )
    return table




def create_condition_select(amx, sample_col):
    columns = np.sort(amx.columns.values)
    columns = np.delete(columns, columns != sample_col)
    select = dcc.Dropdown(
        id = 'condition-select',
        options = [
          {'label' : col, 'value' : col} for col in columns
        ],
        value = columns[0])
    return select



def launch_application(net, gem, amx, vlayers, elayers, sample_col,
    net_name, debug=False):

    """
    Creates the Dash application.

    The Dash application will provide all of the interactive plots, tables and
    filters to interacitvely exploring the network.

    net :  The network dataframe created by the load_network function.

    gem :  The GEM dataframe created by the load_gem function.

    amx : The annotation matrix dataframe created by the load_amx function.

    vlayers : The dataframe containing the 3D coordinates for the nodes.

    elayers : The dataframe containing the 3D coordinates for the edges.

    sample_col : The name of the column in the amx that contains the sample
                 name.

    net_name : The name of the network to display.

    debug:  Set to True to enable Dash debugging support.

    """
    app = dash.Dash()
    app.layout = html.Div([
        # Header Row
        html.Div(className='row', id = "header", children=[
            html.Img(
                src="https://raw.githubusercontent.com/SystemsGenetics/KINC/master/docs/images/kinc.png",
                style={"height" : "55px","display" : "inline-block",
                  "padding" : "0px", "margin" : "0px 10px 0px 10px"
                }),
            html.H1(children="3D Network Explorer",
                style={
                  "display" : "inline-block", "padding" : "10px 0px 0px 0px",
                  "margin" : "0px", "vertical-align" : "top"
                }),
            html.Div(children="Network name: " + net_name,
                style={
                  "padding" : "0px 0px 0px 10px"
                }),
        ]),
        # Graph Row
        html.Div(children=[
            # 3D network plot
            html.Div(children=[
                dcc.Graph(id = 'network-3dview',
                  figure = create_network_plot(net, vlayers, elayers),
                )],
                style = {"width" : "45%", "display" : "inline-block",
                         "border" : "1px solid black", "padding" : "10px",
                         "background-color" : "black", "margin" : "10px"},
            ),
            # 3D Co-Expression scatterplot.
            html.Div(children=[
                html.Div(id='condition-select-box', children=[
                    create_condition_select(amx, sample_col)],
                    style={'padding-bottom' : '10px'}),
                dcc.Graph(id = 'edge-expression-3dview',
                    figure = create_expression_scatterplot(gem, amx, elayers),
                ),
                dcc.Input(
                    id='current-edge', type="number",
                    value=0, style= {'display' : 'none'})],
                style = {"width" : "45%", "display" : "inline-block",
                         "border" : "1px solid black", "padding" : "10px",
                         "margin" : "10px", "vertical-align" : "top"},
            ),
        ]),
        # Table Row
        html.Div(className='row', children=[
            html.Div(id = "edge-table", children=[
                create_dash_edge_table(net),
            ]),
            # html.Div(children = [
            #    html.Pre(id='click-data',
            #        style= {'border': 'thin lightgrey solid',
            #                'overflowX': 'scroll'})],
            # ),
        ]),
    ])

    # Callback when an object in the network plot is clicked.
    @app.callback(
        dash.dependencies.Output('current-edge', 'value'),
        [dash.dependencies.Input('network-3dview', 'clickData')])
    def set_current_edge(clickData):
        if (clickData):
            points = clickData['points']
            found = re.match('^Edge: (.*?) \(co\) (.*?)$', points[0]['text'])
            if (found):
                edge_index = points[0]['customdata']
                return edge_index
        raise dash.exceptions.PreventUpdate

    # Callback to update the co-expression plot when the edge changes.
    @app.callback(
        dash.dependencies.Output('edge-expression-3dview', 'figure'),
        [dash.dependencies.Input('current-edge', 'value'),
         dash.dependencies.Input('condition-select', 'value')])
    def update_expression_plot(current_edge, color_col):
        return create_expression_scatterplot(gem, amx, elayers, color_col, current_edge)

    # Callback for replacing the edge table.
    @app.callback(
         dash.dependencies.Output('edge-table', 'children'),
         [dash.dependencies.Input('current-edge', 'value')])
    def update_edge_table(current_edge):
        return create_dash_edge_table(net, current_edge)

    # Callback for displaying the click data.
    # @app.callback(
    #     dash.dependencies.Output('click-data', 'children'),
    #     [dash.dependencies.Input('network-3dview', 'clickData')])
    # def update_click_data(clickData):
    #     return json.dumps(clickData, indent=2)


    # @app.callback(
    #      dash.dependencies.Output('edge-expression-3dview', 'figure'),
    #      [dash.dependencies.Input('condition-select', 'value'),
    #       dash.dependencies.Input('edge-expression-3dview', 'figure')])
    # def update_condition_select(input, figure):
    #     print(figure)
    #     figure = create_expression_scatterplot(gem, amx, elayers, color_col)

    app.run_server(debug=debug)





def main():
    """
    The main function.

    """

    parser = argparse.ArgumentParser()
    parser.add_argument('--net', dest='net_path', type=str, required=True)
    parser.add_argument('--emx', dest='gem_path', type=str, required=True)
    parser.add_argument('--amx', dest='amx_path', type=str, required=True)
    parser.add_argument('--sample_col', dest='sample_col', type=str, required=False, default='Sample')
    parser.add_argument('--debug', dest='debug', action='store_true', default=False)
    args = parser.parse_args()

    # Make sure the paths exist
    if (not os.path.exists(args.net_path)):
        print ("ERROR: The network file cannot be found: {}".format(args.net_path))
        exit(1)
    if (not os.path.exists(args.gem_path)):
        print ("ERROR: The expression matrix file cannot be found: {}".format(args.gem_path))
        exit(1)
    if (not os.path.exists(args.amx_path)):
        print ("ERROR: The annotation matrix file cannot be found: {}".format(args.amx_path))
        exit(1)

    # Load the input data.
    print("Reading input files...")
    net = load_network(args.net_path)
    gem = load_gem(args.gem_path)
    amx = load_amx(args.amx_path, args.sample_col)


    # Get the filename of the network file minus the extension.
    (net_prefix, net_ext) = os.path.splitext(os.path.basename(args.net_path))


    # Calculate a 2D layout for the network
    print("Calculating 2D layout file:\n  {} \n  This may take awhile if it isn't already present...".format(net_prefix + '.2Dlayout.txt'))
    glayout = calculate_2d_layout(net, net_prefix)

    # Calculate the Z-coorinate positions for the verticies and edges.
    print("Calculating 3D layout...")
    net['Bin'] = bin_edges(net)
    vlayers = get_vertex_zlayers(net, glayout)
    elayers = get_edge_zlayers(net, glayout)

    # Launch the dash application
    print("Launching application...")
    launch_application(net, gem, amx, vlayers, elayers, args.sample_col, net_prefix, args.debug)

    exit(0)



if __name__ == "__main__":
    main()
