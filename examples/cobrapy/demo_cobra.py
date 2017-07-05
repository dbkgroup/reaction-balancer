import cobra
import cobra.test
import requests, json
import cobra_balance as cb

#URL or balancing service
BALANCER_URL = 'http://knime.mib.man.ac.uk:8080/balance'

'''
Demonstration of using the JSON reaction balancing web service from cobrapy
'''

#Load model (see cobrapy tutorials)
model = cobra.test.create_test_model("textbook")
rxn_id = 'PGI'

#For demo purposes we're going to break the reaction stoichiometries
print 'Deliberately disrupting stoichiometries for demo purposes (doubling f6p_c coefficient)'
reaction = model.reactions.get_by_id(rxn_id)
reaction.add_metabolites({'f6p_c':2.0},False)

print 'Original reaction:', reaction.reaction
print ''
model = cb.balance_reaction(model,rxn_id,BALANCER_URL)
new_reaction = model.reactions.get_by_id(rxn_id)
print 'Updated cobra reaction:', new_reaction.reaction

