'''
dbkgroup (c) University of Manchester 2016

reaction-balancer is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
@author:  pauldobson
'''

from flask import Flask, request
from flask_restful import Resource
import balance, json, tempfile
import sbml_service as sbserv
import traceback

#Initialise app
app = Flask(__name__)

#Where the balancing hookup occurs
def list_balancer(listOfRxns):
	results = {}
	for id in listOfRxns:
		rxn = listOfRxns[id]
		try:
			print 'Balancing...'
			is_balanced, was_balanced, balanced_rxn = balance.balance_reac(rxn)
			print '...done'
			print 'Simplifying...'
			balanced_rxn = simplifier(balanced_rxn)
			print '...done'
			if was_balanced:
				msg = 'already balanced'
			else:
				if is_balanced:
					msg = 'brought into balance'
				else:
					msg = 'not brought into balance'
			results[id] = [balanced_rxn,was_balanced,is_balanced,msg]
		except Exception, e:
			#print str(e)
			#traceback.print_exc()
			results[id] = [False, False, False, str(e)]
	return results

#Check minimal reaction stoichiometries are enforced
def simplifier(rxn):
	met_counter = {}
	met_tracker = {}
	#Roll through metabolites adding stoichiometries for equivalents
	for met in rxn:
		formula = met[0]
		charge = met[1]
		name = met[3]
		key = formula + "_" + str(charge) + "_" + name
		stoich = met[2]
		if key in met_counter:
			print 'Naughty! Balancing was required'
			val = met_counter[key]
			met_counter[key] = val + stoich
		else:
			met_counter[key] = stoich
		met_tracker[key] = [formula, charge, False, name]
	#Rebuild reaction
	new_reaction = []
	for met in met_tracker:
		val = met_tracker[met]
		stoich = met_counter[met]
		val[2] = stoich
		new_reaction.append(val)
	return new_reaction

#Endpoint for humans to read
@app.route('/balance', methods=['GET'])
def human():
	message = 	('<html>'
				'<head><title>reaction-balancer</title></head>'
				'<body>'
				'This web service checks and corrects the stoichiometries of chemical reactions.<br><br>'
				'See <a href="https://github.com/metabolicmodelling/reaction-balancer">Github<a/> for full details'
				'</body>'
				'</html>')
	return message

#JSON endpoint
@app.route('/balance/json', methods=['POST'])
def list_json():
	listOfRxns = request.get_json(force=True)
	results = list_balancer(listOfRxns)
	return json.dumps(results)

#SBML endpoint
@app.route('/balance/SBML', methods=['POST'])
def sbml():
	try:
		sbmlString = request.data
	except:
		return 'Could not access SBML data'
	response = sbserv.sbml_balancer(sbmlString)
	return response

if __name__ == '__main__':
	#main()
	app.run(host='0.0.0.0',debug=True, port=8080, threaded=True)