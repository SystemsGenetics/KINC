"""
Visualization script for KINC plain-text co-expression network files.

This script can be  launched as a 'panel' application via:

```
panel serve viz_kinc_network.py --args <file path>
```

Or as a standard python module via:

```
python viz_kinc_network.py <file path>
```

"""

import pandas as pd
import networkx as nx
import holoviews as hv
import param
from textwrap import dedent
import panel as pn
# from methodtools import lru_cache
import argparse
import sys
from bokeh.server.server import Server
from tornado.ioloop import IOLoop
from bokeh.application.handlers import FunctionHandler
from bokeh.application import Application


pn.extension()
hv.extension("bokeh")


class KINCnetworkViz(param.Parameterized):
    """
    Network vizualization application for (plain-text) KINC co-expression files.

    The plain text KINC co-expression file is an edge list.

    """

    data = param.ClassSelector(class_=pd.DataFrame, precedence=-1.0, doc=dedent("""\
    A pandas.DataFrame object sourced from a plain-text KINC co-expression file."""))

    graph = param.Parameter(precedence=-1.0, doc=dedent("""\
    A networkx.Graph object sourced from the `data` parameter. This is populated automatically."""))

    standard_edge_attrs = param.List(default=['Source', 'Target', 'Similarity_Score', 'Interaction', 'Cluster_Index', 'Cluster_Size', 'Samples'],
    precedence=-1.0,
    doc=dedent("""\
    The list of standard edge attributes in the plain-text KINC co-expression files. This should not require
    modification by the user."""))

    conditional_edge_attrs = param.List(precedence=-1.0, doc=dedent("""\
    Conditional edge labels. These are automatically populated during initialization, and should not require
    modification by the user."""))

    node_size = param.Integer(default=5, doc=dedent("""\
    The base size of the nodes."""))

    edge_color = param.Selector(default=None)

    edge_line_width = param.Selector(default=None)

    # bundle_edges = param.Boolean(default=False)

    iterations = param.Integer(default=25)

    graph_opts = param.Dict(default=dict(
        width=800,
        height=800,
        padding=0.1,
        node_line_width=0.1,
        bgcolor="lightgrey",
        edge_cmap="bkr",
        edge_line_width=0.2,
        xaxis=None,
        yaxis=None,
    ), precedence=-1.0)

    def __init__(self, **params):
        super().__init__(**params)
        # Parse the data parameter for conditional edges.
        conditional_edges = self._get_conditional_columns()
        self.conditional_edge_attrs = conditional_edges
        # Populate the selectors.
        self.param["edge_color"].objects = conditional_edges + [None]
        self.param["edge_line_width"].objects = conditional_edges + [None]

        # Build the graph.
        self.graph = nx.from_pandas_edgelist(self.data, source="Source", target="Target", edge_attr=True)

    def _get_conditional_columns(self):
        conditional_columns = [col for col in self.data.columns if col not in self.standard_edge_attrs]
        return conditional_columns

    @classmethod
    def from_filepath(cls, filepath):
        data = pd.read_csv(filepath, sep="\t")
        return cls(data=data)

    # @lru_cache()
    def create_layout(self, frozen_layout_args):
        layout_kwargs = dict(frozen_layout_args)
        return nx.layout.spring_layout(self.graph, **layout_kwargs)

    update_graph = param.Action(lambda self: self.param.trigger('update_graph'),
                                doc="Update view.")

    # TODO: Find a way to enable drag and drop:
    # https://holoviews.org/reference/streams/bokeh/PointDraw.html
    # TODO: Consider enabling bundle_edges
    @param.depends('update_graph', 'edge_color', 'edge_line_width', 'node_size')
    def view(self):
        layout_kwargs = {"iterations": self.iterations,
                         "weight": "Similarity_Score"}
        frozen_layout_kwargs = frozenset(layout_kwargs.items())
        layout = self.create_layout(frozen_layout_kwargs)
        hv_graph = hv.Graph.from_networkx(self.graph, positions=layout)

        opts = {"node_size": self.node_size,
                "edge_color": self.edge_color,
                "edge_line_width": self.edge_line_width,
                **self.graph_opts}

        # if self.bundle_edges:
        #     return bundle_graph(hv_graph.opts(**opts))

        return hv_graph.opts(**opts)

    def panel(self):
        controls = pn.Param(self.param, widgets={
            'update_graph': {'type': pn.widgets.Button, 'button_type': 'primary'},
            'node_size': {'type': pn.widgets.IntSlider, 'start': 1, 'end': 10, 'step': 1, 'value': 5},
        })

        layout = pn.Column("### KINC Network Vizualization",
                           pn.Row(self.view, controls))

        data_table = hv.Table(self.data).opts(width=1200, height=800)

        help_text = generate_help_pane(self)

        return pn.Tabs(("Network", layout),
                       ("Data Table", data_table),
                       ("Documentation", help_text))


def generate_help_pane(self, cutoff=-1.0):
    """
    Generate a panel of the docstrings of the model provided.
    """
    docs = []
    for val in self.param.objects():
        if self.param[val].precedence and self.param[val].precedence >= cutoff:
            docs.append(f"### {self.param[val].label}\n\n" + (self.param[val].doc or ""))
    return pn.pane.Markdown("\n\n".join(docs), width=800)


def main():
    """
    Run the requested visualization. This function is called when this script is called via:
    `python -m kinc_network.py`.

    https://stackoverflow.com/a/43211451
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("filepath")
    args = parser.parse_args()

    instance = KINCnetworkViz.from_filepath(args.filepath).panel()

    def modify_doc(doc):
        doc.add_root(instance.get_root())
        doc.title = "KINC Network Visualization"

    io_loop = IOLoop.current()
    bokeh_app = Application(FunctionHandler(modify_doc))
    server = Server({"/": bokeh_app}, io_loop=io_loop)
    server.start()
    io_loop.add_callback(server.show, "/")
    io_loop.start()


if __name__ == "__main__":
    main()
else:
    KINCnetworkViz.from_filepath(sys.argv[1]).panel().servable()
