import requests, json

#URL or balancing service
BALANCER_URL = 'http://knime.mib.man.ac.uk:8080/balance'

#Function to send a JSON format reaction to the JSON endpoint
def json_balancer(reaction,URL=BALANCER_URL):
	web_service_response = requests.post(URL+'/json', json.dumps(reaction) )
	return web_service_response.json()

#Function to send SBML Level 3 with fbc package to the SBML endpoint
def sbml_balancer(sbmlString,URL=BALANCER_URL):
	web_service_response = requests.post(URL+'/SBML', sbmlString )
	return web_service_response.text

#JSON endpoint
print 'Demonstrating JSON endpoint'
f = open('../data/sample.json', 'r') 				#Load sample JSON
reaction = json.loads(f.read())						#Convert JSON
print ''
print 'Sending:', reaction
print ''
print 'Received:', json_balancer(reaction)			#Hit JSON endpoint and report
print ''

#SBML endpoint
print 'Demonstrating SBML endpoint'
f = open('../data/sample_sbml.xml', 'r') 			#Load sample SBML
sbmlString = f.read()								
print ''
print 'Sending:'
print sbmlString
print ''
print 'Received:'
print sbml_balancer(sbmlString)						#Hit SBML endpoint and report
print ''