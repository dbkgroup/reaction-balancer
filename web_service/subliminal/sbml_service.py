from libsbml import *
import balance, json

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

def serviceFormat_to_sbml(doc,processed,lost):
	print 'serviceFormat_to_sbml'
	model = doc.getModel()
	listOfReactions = model.getListOfReactions()
	for reaction in listOfReactions:
		rid = reaction.getId()
		print 'Reaction ID:', rid
		#Have we processed this id?
		if rid in processed:
			print 'Process', rid
			dense = processed[rid]
			rxn = dense[0]
			was_balanced = dense[1]
			is_balanced = dense[2]
			msg = dense[3]
			if was_balanced:
				note = 'AUTO-BALANCER: already balanced'
			else:
				if is_balanced:
					#Make a mini dictionary of results
					mini = {}
					for c in rxn:
						id = c[3]
						stoich = c[2]
						mini[id] = stoich
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
			print 'Unassessed', rid
			note = 'AUTO-BALANCER: could not be assessed'

		#Update notes
		print 'Set notes'
		if reaction.isSetNotes():
			reaction.appendNotes('<body xmlns="http://www.w3.org/1999/xhtml"><p>' + note + '</p></body>')
		else:
			reaction.setNotes('<body xmlns="http://www.w3.org/1999/xhtml"><p>' + note + '</p></body>')

	print 'Write SBML'
	newSBML = writeSBMLToString(doc)
	print 'Wrote SBML'

	return newSBML

def species_details(model):
	speciesDict = {}
	listOfSpecies = model.getListOfSpecies()
	for species in listOfSpecies:
		id = species.getId()
		splugin = species.getPlugin("fbc")
		status = True
		#Get charge
		try:
			charge = splugin.getCharge()
		except:
			status = False
			charge = False
		#Get formula
		try:
			chemicalFormula = splugin.getChemicalFormula()
		except:
			status = False
			chemicalFormula = False
		#print id, [chemicalFormula, charge, status]
		speciesDict[id] = [chemicalFormula, charge, status]
	return speciesDict

def sbml_to_serviceFormat(model):
	print 'sbml_to_service'
	#Get species
	speciesDict = species_details(model)

	#Get reactions
	processable_rxns = {}		#Can be sent to balancer
	lost_rxns = {}				#Cannot be sent to balancer (some components underspecified)
	listOfReactions = model.getListOfReactions()
	for rxn in listOfReactions:
		dense = []
		reactionStatus = True
		rid = rxn.getId()
		reactants = rxn.getListOfReactants() #list of speciesReferences
		for r in reactants:
			#Species id
			mid = r.getSpecies()
			#Pull species details
			details = speciesDict[mid]
			formula = details[0]
			charge = details[1]
			status = details[2]
			#Get stoichiometry
			stoich = r.getStoichiometry()
			if status:
				row = [formula,charge,-1*stoich,mid]
				dense.append(row)
			else:
				reactionStatus = False
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
		if reactionStatus:
			processable_rxns[rid] = dense
		else:
			lost_rxns[rid] = False
	return processable_rxns, lost_rxns

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
			print len(processable), 'reactions can be processed'
			print len(lost), 'reactions cannot be processed'
			print 'Ran sbml_to_service'
		except:
			return 'Could not convert SBML'

		#Hit balancer
		try:
			processed = list_balancer(processable)
			print 'Hit list_balancer'
		except:
			return 'Problem hitting balancer'

		#Push results into model
		try:
			new_sbml = serviceFormat_to_sbml(doc,processed,lost)
			print 'Ran serviceFormat_to_sbml'
			response = new_sbml
			print 'Processed:', processed
			return json.dumps(processed)
		except Exception as e:
			return 'Could not update SBML file with balance results', str(e)

	return False
