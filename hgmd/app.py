import os

from flask import (
    Flask, flash, request, redirect, url_for, render_template,
    send_from_directory
)
from werkzeug.utils import secure_filename

from . import hgmd as md

UPLOAD_FOLDER = 'data/input/'
OUTPUT_FOLDER = 'data/output/'
ALLOWED_EXTENSIONS = set(['txt'])

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['OUTPUT_FOLDER'] = OUTPUT_FOLDER


@app.route('/', methods=['GET', 'POST'])
def hello_world():
    if request.method == 'POST':
        if 'file' not in request.files:
            flash("No file part")
            return redirect(request.url)
        file = request.files['file']
        if file.filename == '':
            flash("No selected file")
            return redirect(request.url)
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            return redirect(url_for('uploaded_file', filename=filename))
    return render_template('app.html')


@app.route('/process', methods=['GET', 'POST'])
def process():
    print("Reading data...")
    cell_data = md.get_cell_data(
        marker_path=(UPLOAD_FOLDER + 'markers.txt'),
        tsne_path=(UPLOAD_FOLDER + 'tsne.txt'),
        cluster_path=(UPLOAD_FOLDER + 'cluster.txt')
    )
    X = 3
    L = 15
    min_exp_ratio = 0.4
    plot_pages = 10
    plot_genes = 10

    # Enumerate clusters and process each individually in its own folder.
    # pair_data also contains singleton data, but singleton is just
    # singleton.
    clusters = cell_data['cluster'].unique()
    clusters.sort()
    for cluster in range(1, 2):
        print("Processing cluster " + str(cluster) + "...")
        cluster_path = OUTPUT_FOLDER + "/cluster_" + str(cluster) + "/"
        os.makedirs(cluster_path, exist_ok=True)
        print("Testing singletons...")
        singleton_data = md.singleton_test(cell_data, cluster, X, L)
        print("Testing pairs...")
        pair_data = md.pair_test(
            cell_data, singleton_data, cluster, min_exp_ratio
        )
        print("Calculating true positive/negative rates...")
        singleton_data, pair_data = md.find_TP_TN(
            cell_data, singleton_data, pair_data, cluster
        )
        print("Saving to CSV...")
        singleton_data.to_csv(cluster_path + "singleton_data.csv")
        pair_data.to_csv(cluster_path + "pair_data.csv")
        print("Done.")
        print("Plotting true positive/negative rates...")
        md.make_TP_TN_plots(
            cell_data, singleton_data, pair_data, plot_genes,
            pair_path=(cluster_path + "TP_TN_plot.pdf"),
            singleton_path=(cluster_path + "singleton_TP_TN_plot.pdf")
        )
        print("Done.")
        print("Plotting discrete expression...")
        md.make_discrete_plots(
            cell_data, singleton_data, pair_data, plot_pages,
            path=(cluster_path + "discrete_plots.pdf"),
        )
        print("Done.")
        print("Plotting continuous expresssion...")
        md.make_combined_plots(
            cell_data, singleton_data, pair_data, plot_genes,
            pair_path=(cluster_path + "combined_plot.pdf"),
            singleton_path=(cluster_path + "singleton_combined_plot.pdf")
        )
        print("Done.")

    print("All set!!! Enjoy your PDFs and CSVs.")
    return redirect(url_for('output_file', filename='cluster_1/pair_data.csv'))


@app.route('/uploads/<filename>')
def uploaded_file(filename):
    return send_from_directory(app.config['UPLOAD_FOLDER'], filename)


@app.route('/output/<filename>')
def output_file(filename):
    return send_from_directory(app.config['OUTPUT_FOLDER'], filename)


def allowed_file(filename):
    return '.' in filename and \
        filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS
