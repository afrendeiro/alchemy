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
from bokeh.templates import RESOURCES
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
    from collections import OrderedDict

    if request.method == 'POST':
        render_template('pages/process.html')

        input_file = request.files['input_file']
        background_file = request.files['background_file']
        if input_file and allowed_file(input_file.filename):
            input_filename = secure_filename(input_file.filename)
            input_file.save(os.path.join(app.config['UPLOAD_FOLDER'], input_filename))

            if background_file and allowed_file(background_file.filename):
                background_filename = secure_filename(background_file.filename)
                background_file.save(os.path.join(app.config['UPLOAD_FOLDER'], background_filename))

            # Assemble arguments
            # input_file
            args = [
                os.path.join(app.config['UPLOAD_FOLDER'], input_filename)]
            # background
            if background_file and allowed_file(background_file.filename):
                args += [
                    "-b",
                    os.path.join(app.config['UPLOAD_FOLDER'], background_filename)]
            # input_type
            args += [
                "-t",
                request.form.get('input_type')]
            # delimiter
            args += [
                "-d",
                request.form.get('delimiter')]
            # header
            # annotation_output
            # background_annotation_output
            # enrichment_output
            # enrich
            # processors
            # tmp_dir

            # start processing data
            annotated_query, enriched_query = alchemy.main(args)

            # return results

            # Bokeh graph
            fig = figure(title="Volcano plot")
            fig.scatter(enriched_query.oddsratio, enriched_query.corrected_pvalue, color="#000000")
            fig.logo = None
            fig.toolbar_location = None

            # Configure resources to include BokehJS inline in the document.
            plot_resources = RESOURCES.render(
                js_raw=INLINE.js_raw,
                css_raw=INLINE.css_raw,
                js_files=INLINE.js_files,
                css_files=INLINE.css_files,
            )
            script, div = components(fig)
            html = flask.render_template(
                'pages/results.html',
                plot_script=script,
                plot_div=div,
                plot_resources=plot_resources,
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
