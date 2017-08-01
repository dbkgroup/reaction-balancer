echo 'Using cURL with JSON'
curl -H "Content-Type: application/json" -X POST -d @../data/sample.json http://knime.mib.man.ac.uk:8080/balance/json

echo 'Using cURL with SBML'
curl -H "Content-Type: application/xml" -X POST -d @../data/sample_sbml.xml http://knime.mib.man.ac.uk:8080/balance/SBML
