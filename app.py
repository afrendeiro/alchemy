#!/usr/bin/env python

"""
This is the Flask app that runs the GUI to the alchemy package.
"""

import os

# Flask
import flask
from flask import Flask, request, render_template
from werkzeug import secure_filename

# Bokeh
from bokeh.embed import components
from bokeh.plotting import figure
from bokeh.resources import INLINE
from bokeh.util.string import encode_utf8

# Flask settings/configs
UPLOAD_FOLDER = 'uploads'
ALLOWED_EXTENSIONS = set(['txt', 'csv'])

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS


@app.route("/", methods=['GET'])
def home():
    if request.method == 'GET':
        return render_template('pages/home.html', title="home")


@app.route("/about", methods=['GET'])
def about():
    if request.method == 'GET':
        return render_template('pages/about.html', title="About")


@app.route("/register", methods=['GET'])
def register():
    if request.method == 'GET':
        return render_template('forms/register.html', title="Register")


@app.route("/login", methods=['GET'])
def login():
    if request.method == 'GET':
        return render_template('forms/login.html', title="Login")


@app.route("/forgot", methods=['GET'])
def forgot():
    if request.method == 'GET':
        return render_template('forms/forgot.html', title="Login")


@app.route("/upload", methods=['GET'])
def upload():
    if request.method == 'GET':
        return render_template('pages/upload.html', files=os.listdir(app.config['UPLOAD_FOLDER'],))


@app.route("/process", methods=['GET', 'POST'])
def process():
    """
    Receive the upload smiles, process and retrieve results.
    """
    from alchemy import alchemy
    # filename = "t.csv"

    if request.method == 'POST':
        render_template('pages/process.html')
        file = request.files['file']
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))

            # start processing data
            (annotated_query, enriched_query) = alchemy.main([
                os.path.join(app.config['UPLOAD_FOLDER'], filename),  # input_file
                # background
                # input_type
                # delimiter
                # header
                # annotation_output
                # background_annotation_output
                # enrichment_output
                # enrich
                # processors
                # tmp_dir
            ])

            # return results
            from collections import OrderedDict

            # Bokeh graph
            # Create a polynomial line graph
            fig = figure(title="Polynomial")
            fig.line(range(5), range(5), color="#000000", line_width=2)
            fig.logo = None
            fig.toolbar_location = None

            # Configure resources to include BokehJS inline in the document.
            # For more details see:
            #   http://bokeh.pydata.org/en/latest/docs/reference/resources_embedding.html#bokeh-embed
            js_resources = INLINE.render_js()
            css_resources = INLINE.render_css()

            # For more details see:
            #   http://bokeh.pydata.org/en/latest/docs/user_guide/embedding.html#components
            script, div = components(fig, INLINE)
            html = flask.render_template(
                'pages/results.html',
                plot_script=script,
                plot_div=div,
                js_resources=js_resources,
                css_resources=css_resources,
                compounds=OrderedDict(zip(
                    enriched_query.index, enriched_query.corrected_pvalue))
            )

            return encode_utf8(html)

    return render_template('pages/home.html')


# Default port:
if __name__ == '__main__':
    app.secret_key = 'super secret key'
    app.config['SESSION_TYPE'] = 'filesystem'

    app.debug = True
    app.run()
