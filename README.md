# Automatic chemical reaction balancer
A web service for automatically balancing chemical reaction stoichiometries.

A live web service is available at ...

You will need to install Docker to run the service locally.

## How to set up the service from Github
Clone or download the repository. From inside the repository...

Run: ```docker build -t balancer .```

Run: ```docker run --rm -it -p XXXX:8080 balancer``` (where XXXX is the port you wish to expose on).

## How to access the service
Two endpoints are defined for JSON and SBML. Both use POST.

### JSON
The JSON endpoint at ```http://<URL>:<PORT>/balance/json``` accepts a reaction in the following format...
```
{ <reaction id> : [	{"formula":<formula>, "charge":<charge>, "stoichiometry":<stoichiometry>, "name":<molecularSpecies_name>]  ]}
```
An example follows (incomplete gross photosynthetic reaction)...
```
{ "photosynthesis" :
  [	{"formula":"CO2", "charge":0, "stoichiometry":-1.0, "name":"carbon dioxide"},
    {"formula":"C6H12O6", "charge":0, "stoichiometry":1.0, "name":"glucose"},
    {"formula":"O2", "charge":0, "stiochiometry":1.0, "name":"oxygen"}
  ]}
```
The response is as follows...
```{ <reaction_id> : [ [reaction], False, True, "brought into balance"  ]  }```
...where the list elements contain the reaction, whether or not it was originally balanced, whether or not it is now balanced, and a message from the balancer.

For example...
```
{ "photosynthesis" :
  {"reaction":
    [ {"formula":"H2O", "charge":0, "stoichiometry":-6.0, "name":"water"},
      {"formula":"CO2", "charge":0, "stoichiometry":-6.0, "name":"carbon dioxide"},
      {"formula":"C6H12O6", "charge":0, "stoichiometry":1.0, "name":"glucose"},
      {"formula":"O2", "charge":0, "stoichiometry":6.0, "name":"oxygen"}
    ],
   "was_balanced":False,
   "is_balanced":True,
   "message":"brought into balance"
  ]
}

```

### SBML
The SBML endpoint at ```http://<URL>:<PORT>/balance/SBML``` expects SBML Level 3 with the 'fbc' package used to store molecular species' ```chemicalFormula``` and ```charge```. Earlier SBML levels will be rejected by the service. Use POST to send SBML as data.

Existing balanced reactions maintain original stoichiometries and tagged (in ```Notes```) as 'already balanced'.
Unbalanced reactions brought into balance are tagged as 'brought into balance'.
Unbalanced reactions not brought into balance are tagged as 'not brought into balance'.
Reactions that could not be processed because one or more components lacked formula or charge are tagged 'could not be assessed'.
