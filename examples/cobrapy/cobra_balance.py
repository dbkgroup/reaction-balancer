import cobra
import cobra.test
import requests, json


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
def json_balancer(reaction,URL):
	web_service_response = requests.post(URL+'/json', json.dumps(reaction) )
	return web_service_response.json()

#Wrap entire conversion
def balance_reaction(model,rxn_id,URL):
	print 'Converting', rxn_id, 'to web service JSON'
	service_format = cobra_to_json(rxn_id,model)

	print 'Using reaction balancing web service'
	balanced = json_balancer(service_format,URL)

	print 'Updating cobra model'
	model = json_to_cobra(balanced,model)

	return model



