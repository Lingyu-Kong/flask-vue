from flask import Flask, render_template, jsonify, request
from flask_uploads import UploadSet, configure_uploads, DATA
from random import *
from flask_cors import CORS
import os
from ase.io import read, write
from ase.calculators.emt import EMT

app = Flask(__name__,
            static_folder = "./dist/static",
            template_folder = "./dist")
cors = CORS(app, resources={r"/api/*": {"origins": "*"}})

custom_extensions = ('xyz', 'res', 'vasp', '.cif')  
files = UploadSet('files', extensions=custom_extensions) 
app.config['UPLOADED_FILES_DEST'] = 'uploads'  
configure_uploads(app, files)

@app.route('/api/upload', methods=['POST'])  
def upload():  
    if 'file' not in request.files:  
        return jsonify({'error': 'No file uploaded'}), 400  
  
    uploaded_file = request.files['file']  
    filename = files.save(uploaded_file)  
  
    # 读取文件内容并进行计算  
    file_path = os.path.join(app.config['UPLOADED_FILES_DEST'], filename)  
    result = EMT_calc(file_path)  
  
    # 删除临时文件  
    os.remove(file_path)  
  
    # 返回计算结果  
    return jsonify({
        "energy": result["energy"],
        "forces": result["forces"].tolist(),
        "stress": result["stress"].tolist(),
    })  
  
def EMT_calc(file_path):  
    # 根据文件内容进行计算，并返回结果  
    # 这里只是一个示例，您需要根据实际需求实现具体的计算逻辑  
    atoms = read(file_path)
    print("Start EMT calculation")
    calc = EMT()
    atoms.set_calculator(calc)
    results = {
        "energy": atoms.get_potential_energy(),
        "forces": atoms.get_forces(),
        "stress": atoms.get_stress()
    }
    return results
    

@app.route('/api/random')
def random_number():
    response = {
        'randomNumber': randint(1, 100)
    }
    return jsonify(response)

@app.route('/', defaults={'path': ''})
@app.route('/<path:path>')
def catch_all(path):
    return render_template("index.html")