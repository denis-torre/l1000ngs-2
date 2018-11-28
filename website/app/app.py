#################################################################
#################################################################
############### Dubois RNA-seq Analysis #########################
#################################################################
#################################################################
##### Author: Denis Torre
##### Affiliation: Ma'ayan Laboratory,
##### Icahn School of Medicine at Mount Sinai

#######################################################
#######################################################
########## 1. App Configuration
#######################################################
#######################################################

#############################################
########## 1. Load libraries
#############################################
##### 1. Flask modules #####
from flask import Flask, render_template, url_for, request
import os, json
import pandas as pd
entry_point = '/L1000NGS'
app = Flask(__name__, static_url_path='/app/static')

##### 2. Prefix middleware #####
class PrefixMiddleware(object):

	def __init__(self, app, prefix=''):
		self.app = app
		self.prefix = prefix

	def __call__(self, environ, start_response):
		if environ['PATH_INFO'].startswith(self.prefix):
			environ['PATH_INFO'] = environ['PATH_INFO'][len(self.prefix):]
			environ['SCRIPT_NAME'] = self.prefix
			return self.app(environ, start_response)
		else:
			start_response('404', [('Content-Type', 'text/plain')])
			return ["This url does not belong to the app.".encode()]
app.wsgi_app = PrefixMiddleware(app.wsgi_app, prefix=entry_point)

#######################################################
#######################################################
########## 2. Routes
#######################################################
#######################################################
##### Handles routes used to generate notebooks.

##################################################
########## 2.1 Webpages
##################################################

#############################################
########## 1. Home
#############################################
### Landing page for the website. Links to analyze() endpoint.
### Links to: analyze().

@app.route('/')
def index():

	# Read samples
	return render_template('index.html')

#############################################
########## 2. Metadata API
#############################################

@app.route('/api/metadata')
def metadata_api():
	
	# Read data
	dataframe = pd.read_table('app/static/data/signature_metadata.txt').head(50)
	dataframe['checkbox'] = ''

	# Return JSON
	return json.dumps({'data': dataframe.to_dict(orient='records')})

#############################################
########## 3. Results
#############################################

@app.route('/analyze', methods=['GET', 'POST'])
def analyze():

	# Get RID list
	# rid_list = request.form.get('rid_list').split(',')
	rid_list = ['AML001_CD34_6H:BRD-K43389675:10', 'AML001_PC3_6H:BRD-A19037878:0.37037', 'AML001_PC3_6H:BRD-A19037878:1.11111']

	# Return JSON
	return str(rid_list)
