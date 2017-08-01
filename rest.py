'''
NaCTeM/DBK Group (c) University of Manchester 2017

reaction-balancer is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
@author:  pauldobson
'''

from flask import Flask, request
from flask_restful import Resource, Api
from flask_cors import CORS, cross_origin
from libsbml import *
#from subliminal import balance
import json

#Initialise app
app = Flask(__name__, static_url_path='')
CORS(app)
api = Api(app)

#----------------------------------------------------------
#Where the balancing hookup occurs
def list_balancer(listOfRxns):
	results = {}
	for id in listOfRxns:
		rxn = listOfRxns[id]
		try:
			is_balanced, was_balanced, balanced_rxn = balance.balance_reac(rxn,[('H', 1, 'h_balancer'), ('H2O', 0, 'h2o_balancer')])
			if was_balanced:
				msg = 'already balanced'
			else:
				if is_balanced:
					msg = 'brought into balance'
				else:
					msg = 'not brought into balance'
			results[id] = {"reaction":balanced_rxn,"was_balanced":was_balanced,"is_balanced":is_balanced,"message":msg}
			print 'List balancer: okay'
		except Exception, e:
			results[id] = {"reaction":rxn, "was_balanced":False, "is_balanced":False, "message":False}
			print 'List balancer: error', str(e)
	return results
#----------------------------------------------------------

#----------------------------------------------------------
#SBML functions

#Extract species' formulae and charges from SBML
def species_details(model):
	speciesDict = {}
	listOfSpecies = model.getListOfSpecies()
	for species in listOfSpecies:
		id = species.getId()
		splugin = species.getPlugin("fbc")
		status = True	#Set flag for reaction specification (if this flag is not switched all components were sufficiently specified)
		#Get charge
		try:
			charge = splugin.getCharge()
		except:
			status = False	#Switch reaction specification flag
			charge = False
		#Get formula
		try:
			chemicalFormula = splugin.getChemicalFormula()
		except:
			status = False	#Switch reaction specification flag
			chemicalFormula = False
		#print id, [chemicalFormula, charge, status]
		speciesDict[id] = [chemicalFormula, charge, status]
	return speciesDict

#Convert SBML to service format
def sbml_to_serviceFormat(model):
	print 'sbml_to_service'

	#Get species information from SBML
	speciesDict = species_details(model)

	#Split reactions into processable (all components sufficiently specified) and lost (at least one component insufficiently specified)
	#sufficiently specified means both formula and charge are reported using the fbc package
	processable_rxns = {}		#Can be sent to balancer
	lost_rxns = {}				#Cannot be sent to balancer (some components underspecified)

	#Get reactions from listOfReactions	
	listOfReactions = model.getListOfReactions()
	for rxn in listOfReactions:
		dense = []
		reactionStatus = True
		rid = rxn.getId()	#Reaction identifier
		
		#Get reactants from listOfReactants
		reactants = rxn.getListOfReactants() 	#list of speciesReferences
		for r in reactants:
			#Species id
			mid = r.getSpecies()
			#Pull species details
			details = speciesDict[mid]
			formula = details[0]
			charge = details[1]
			status = details[2]	#Status flags whether or not a molecular species is sufficiently specified for the balancer
			#Get stoichiometry
			stoich = r.getStoichiometry()
			if status:	#If sufficiently specified add to list
				row = [formula,charge,-1*stoich,mid]
				dense.append(row)
			else:		#Else bail out and flag reaction as lost (cannot be balanced)
				reactionStatus = False
				
		#Get products from listOfProducts
		products = rxn.getListOfProducts() 
		for p in products:
			#Species id
			mid = p.getSpecies()
			#Pull species details
			details = speciesDict[mid]
			formula = details[0]
			charge = details[1]
			status = details[2]
			#Get stoichiometry
			stoich = p.getStoichiometry()
			if status:
				row = [formula,charge,stoich,mid]
				dense.append(row)
			else:
				reactionStatus = False

		#Validate reaction
		if reactionStatus: 	#If all components are sufficiently specified add reaction to the processable list
			processable_rxns[rid] = dense
		else:				#If some components are insufficiently specified add reaction to the list list
			lost_rxns[rid] = False
			
	return processable_rxns, lost_rxns

#Convert service format to SBML
def serviceFormat_to_sbml(doc,processed,lost):
	print 'serviceFormat_to_sbml'
	model = doc.getModel()

	#Add balance compartment
	model.createCompartment()
	c1 = model.createCompartment()
	c1.setId('BALANCER')
	c1.setConstant(True)
	c1.setSize(1)
	c1.setSpatialDimensions(3)
	c1.setUnits('litre')
	
	#Add balance metabolites
	s1 = model.createSpecies()
	s1.setId('h2o_balancer')
	s1.setCompartment('BALANCER')
	s1.setConstant(False)
	s1.setInitialAmount(5)
	s1.setSubstanceUnits('mole')
	s1.setBoundaryCondition(True)
	s1.setHasOnlySubstanceUnits(False)
	splugin1 = s1.getPlugin("fbc")
	splugin1.setCharge(0)
	splugin1.setChemicalFormula("H2O")

	s2 = model.createSpecies()
	s2.setId('h_balancer')
	s2.setCompartment('BALANCER')
	s2.setConstant(False)
	s2.setInitialAmount(5)
	s2.setSubstanceUnits('mole')
	s2.setBoundaryCondition(True)
	s2.setHasOnlySubstanceUnits(False)
	splugin2 = s2.getPlugin("fbc")
	splugin2.setCharge(1)
	splugin2.setChemicalFormula("H")
	
	#Process reactions
	listOfReactions = model.getListOfReactions()
	for reaction in listOfReactions:
		rid = reaction.getId()
		#print 'Reaction ID:', rid
		#Have we processed this id?
		if rid in processed:
			#print 'Process', rid
			dense = processed[rid]
			rxn = dense[0]
			was_balanced = dense[1]
			is_balanced = dense[2]
			msg = dense[3]
			if was_balanced:
				note = 'AUTO-BALANCER: already balanced'
			else:
				if is_balanced:
					#Make a mini dictionary of results and add in balancer metabolites
					mini = {}
					for c in rxn:
						id = c[3]
						stoich = c[2]
						mini[id] = stoich
						if id is 'h2o_balancer':
							if stoich < 0:
								reactant = reaction.createReactant()
								reactant.setSpecies("h2o_balancer")
								reactant.setStoichiometry(-1*stoich)
								reactant.setConstant(False)
							else:
								product = reaction.createProduct()
								product.setSpecies("h2o_balancer")
								product.setStoichiometry(stoich)
								product.setConstant(False)
						if id is 'h_balancer':
							if stoich < 0:
								reactant = reaction.createReactant()
								reactant.setSpecies("h_balancer")
								reactant.setStoichiometry(-1*stoich)
								reactant.setConstant(False)
							else:
								product = reaction.createProduct()
								product.setSpecies("h_balancer")
								product.setStoichiometry(stoich)
								product.setConstant(False)
						
					#Push results into model
					reactants = reaction.getListOfReactants()
					for reactant in reactants:
						mid = reactant.getSpecies()
						new_stoich = mini[mid]
						reactant.setStoichiometry(new_stoich)
					products = reaction.getListOfProducts()
					for product in products:
						mid = product.getSpecies()
						new_stoich = mini[mid]
						product.setStoichiometry(new_stoich)	
					note = 'AUTO-BALANCER: brought into balance'
				else:
					note = 'AUTO-BALANCER: not brought into balance'
		#Report that the balancer could not be applied
		else:
			#print 'Unassessed', rid
			note = 'AUTO-BALANCER: could not be assessed'

		#Update notes
		#print 'Set notes'
		if reaction.isSetNotes():
			reaction.appendNotes('<body xmlns="http://www.w3.org/1999/xhtml"><p>' + note + '</p></body>')
		else:
			reaction.setNotes('<body xmlns="http://www.w3.org/1999/xhtml"><p>' + note + '</p></body>')

	#print 'Write SBML'
	newSBML = writeSBMLToString(doc)
	#print 'Wrote SBML'

	return newSBML

#Run SBML balancer
def sbml_balancer(string):
	#Read SBML
	reader = SBMLReader()
	try:
		doc = reader.readSBMLFromString(string)
		print 'Successful readSBMLFromString'
	except:
		return 'Could not readSBMLFromString'

	#Check for errors
	errors = doc.getNumErrors();
	if doc.getNumErrors() > 0:
		if doc.getError(0).getErrorId() == XMLFileUnreadable:
			print 'Unreadable XML'
			return 'Unreadable XML'
		elif doc.getError(0).getErrorId() == XMLFileOperationError:
			print 'Operation error on XML'
			return 'Operation error on XML'
		else:
			print 'General cockup on the SBML reading front'
			return 'General cockup on the SBML reading front'

	#Get model
	model = doc.getModel()
	if model == None:
		print 'Unable to create model object'
		return 'Unable to create model object'
	
	#Check for fbc plugin
	mplugin = model.getPlugin("fbc")
	if mplugin == None:
		print 'Plugin error: the SBML service expects formulae and charges to be represented using the fbc package'
		return 'Plugin error: the SBML service expects formulae and charges to be represented using the fbc package'
	else:
		#Turn the SBML into service format
		try:
			processable, lost = sbml_to_serviceFormat(model)
			#print processable
			print len(processable), 'reactions can be processed'
			print len(lost), 'reactions cannot be processed'
			print 'Ran sbml_to_service'
		except:
			return 'Could not convert SBML'

		#Hit balancer
		try:
			processed = list_balancer(processable)
			#print processed
			print 'Hit list_balancer'
		except:
			return 'Problem hitting balancer'

		#Push results into model
		try:
			new_sbml = serviceFormat_to_sbml(doc,processed,lost)
			print 'Ran serviceFormat_to_sbml'
			response = new_sbml
			#print 'Processed:', processed
			return new_sbml
		except Exception as e:
			return 'Could not update SBML file with balance results', str(e)

	return False
#----------------------------------------------------------

#----------------------------------------------------------
#ENDPOINTS

#Landing page for human interface
@app.route('/')
def index():
	return app.send_static_file('index.html')

#Endpoint to be used by human interface
@app.route('/balance/manual', methods=['POST'])
def manual():
	print 'Human UI'
	content = request.get_json(force=True)
	tidyList = []
	for c in content:
		if str(c[0]) is not '':
			tidyComponent = [ str(c[0]), int(c[1]), float(c[2]), str(c[3]) ]
			tidyList.append(tidyComponent)

	#Use the user-supplied reaction
	if len(tidyList) > 0:
		serviceFormat = { "human" : tidyList}
	#Use the default reaction
	else:
		print 'Running default...'
		serviceFormat = { "human" : [["CO2",0,-1.0,"carbon dioxide"],["C6H12O6",0,1.0,"glucose"],["O2",0,1.0,"oxygen"]] }

	#Run reaction
	print 'Hitting balancer...'
	results = list_balancer(serviceFormat)
	print '...complete'
	
	#Process results
	print 'Converting result...'
	response = results["human"]
	reaction = response["reaction"]
	was_balanced = response["was_balanced"]
	is_balanced = response["is_balanced"]
	msg = response["message"]
	print reaction
	
	#Make reaction into string format
	reactants = []
	products = []
	for component in reaction:
		print component
		formula = component[0]
		charge = component[1]
		stoich = int(component[2])
		name = component[3]
		if stoich < 0:
			stoich = stoich * -1
			elem = str(stoich) + '.' + name
			reactants.append(elem)
		else:
			elem = str(stoich) + '.' + name
			products.append(elem)
	reactantString = ' + '.join(reactants)
	productString = ' + '.join(products)
	reactionString = reactantString + '  ->  ' + productString
	print reactionString
	
	return reactionString

#JSON endpoint
@app.route('/balance/json', methods=['POST'])
def list_json():
	print 'JSON balancer'
	verboseListOfRxns = request.get_json(force=True)
	print 'Incoming:', verboseListOfRxns

	#Convert verbose JSON to Subliminal format
	print 'Convert verbose JSON to Subliminal format'
	listOfRxns = {}
	for rid in verboseListOfRxns:
		print rid, 
		rxn = verboseListOfRxns[rid]
		print rxn
		
		#Set of reaction components in lite subliminal format
		subliminal = []
		for mol in rxn:
			formula = str(mol["formula"])
			charge = str(mol["charge"])
			stoichiometry = mol["stoichiometry"]
			name = str(mol["name"])
			row = [formula,charge,stoichiometry,name]
			subliminal.append(row)
		
		#Add set to listOfRxns
		listOfRxns[rid] = subliminal
	print 'Subliminal (all):', listOfRxns

	#Run balance
	print 'Run balancer'
	newListOfRxns = list_balancer(listOfRxns)
	print 'Results:', newListOfRxns

	#Convert Subliminal format to more verbose JSON
	print 'Convert Subliminal format to verbose JSON'
	results = {}
	for rid in newListOfRxns:
		print rid
		result = newListOfRxns[rid]
		was_balanced = result["was_balanced"]
		is_balanced = result["is_balanced"]
		message = result["message"]
		reaction = result["reaction"]

		#Get the list of reaction components
		verbose_reaction = []
		for mol in reaction:
			formula = mol[0]
			charge = mol[1]
			stoichiometry = mol[2]
			name = mol[3]
			verbose_reaction.append( {"formula":formula,"charge":charge,"stoichiometry":stoichiometry,"name":name} )

		#Update results list
		results[rid] = {"reaction":verbose_reaction,"was_balanced":was_balanced,"is_balanced":is_balanced,"message":message}

	return json.dumps(results)

#SBML endpoint
@api.representation('application/xml')
@app.route('/balance/SBML', methods=['POST'])
def sbml():
	try:
		sbmlString = request.data
	except:
		return 'Could not access SBML data'
	response = sbml_balancer(sbmlString)
	return response
#----------------------------------------------------------

if __name__ == '__main__':
	#main()
	app.run(host='0.0.0.0',debug=False, port=8080, threaded=True)