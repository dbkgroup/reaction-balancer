import cobra
import cobra.test
import requests, json

#URL or balancing service
BALANCER_URL = 'http://knime.mib.man.ac.uk:8080/balance'

'''
Demonstration of using the JSON reaction balancing web service from cobrapy
'''

#Convert a cobrapy reaction to web service format
def cobra_to_json(rxn_id,model):
	try:
		#Container for the web service
		service_json = []

		#Get reaction from model
		rxn = model.reactions.get_by_id(rxn_id)

		#Get reactants
		reactants = rxn.reactants
		for r in reactants:
			#Get metabolite properties
			metabolite = model.metabolites.get_by_id(str(r))
			stoich = rxn.get_coefficient(str(r))
			service_json.append( [metabolite.formula, metabolite.charge, stoich, str(r)] )

		#Get products
		products = rxn.products
		for p in products:
			#Get metabolite properties
			metabolite = model.metabolites.get_by_id(str(p))
			stoich = rxn.get_coefficient(str(p))
			service_json.append( [metabolite.formula, metabolite.charge, stoich, str(p)] )

		#Assemble JSON
		service_format = {rxn_id:service_json}
		return service_format
	except:
		return 'Error converting reaction to web service format (do all reaction components have formula and charge?)'

#Push web service results back into the cobra model
def json_to_cobra(service_response,model):
	print ''
	#Service response contains a list of reactions. Process one at a time.
	for rxn_id in service_response:

		#Get original reaction
		old_cobra_reaction = model.reactions.get_by_id(rxn_id)
		print 'Original cobra reaction:', old_cobra_reaction.reaction

		#Process JSON response
		new_reaction, was_balanced, is_balanced, msg = service_response[rxn_id]
		if not was_balanced:
			if is_balanced:
				reaction = model.reactions.get_by_id(rxn_id)
				for met in new_reaction:
					formula, charge, stoich, met_id = met
					stoich = reaction.add_metabolites({str(met_id):stoich},False)
	return model

#Function to send a JSON format reaction to the JSON endpoint
def json_balancer(reaction,URL=BALANCER_URL):
	web_service_response = requests.post(URL+'/json', json.dumps(reaction) )
	return web_service_response.json()


#Load model (see cobrapy tutorials)
model = cobra.test.create_test_model("textbook")

#For demo purposes we're going to break the reaction stoichiometries
print 'Deliberately disrupting stoichiometries for demo purposes (doubling f6p_c coefficient)'
reaction = model.reactions.get_by_id('PGI')
reaction.add_metabolites({'f6p_c':2.0},False)
print 'Original reaction:', reaction.reaction
print ''

#Convert cobrapy reaction 'PGI' to JSON for web service
service_format = cobra_to_json("PGI",model)
print 'Converting PGI to web service JSON'

print 'Using reaction balancing web service'
balanced = json_balancer(service_format)

print 'Updating cobra model'
model = json_to_cobra(balanced,model)

new_reaction = model.reactions.get_by_id('PGI')
print 'Updated cobra reaction:', new_reaction.reaction