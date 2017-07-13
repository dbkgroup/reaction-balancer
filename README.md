# Automatic chemical reaction balancer
A web service for automatically balancing chemical reaction stoichiometries. You will need to install Docker to run the service locally.

## From Docker
```docker run -p XXXX:8080 compsysbio/reaction-balancer``` (where XXXX is the port you wish to expose on).

## How to set up the service from Github
Clone or download the repository.

Move into the web_service folder.

Run: ```docker build -t balancer .```

Run: ```docker run --rm -it -p XXXX:8080 balancer``` (where XXXX is the port you wish to expose on).

## How to access the service
Two endpoints are defined for JSON and SBML. Both use POST.

### JSON
The JSON endpoint at ```http://<URL>:<PORT>/balance/json``` accepts a reaction in the following format...
```
{ <reaction id> : [	[formula, charge, stoichiometry, <molecularSpecies_id>]  ]}
```
An example follows (incomplete gross photosynthetic reaction)...
```
{ "photosynthesis" :
  [	["CO2", 0, -1.0, "carbon dioxide"],
    ["C6H12O6", 0, 1.0, "glucose"],
    ["O2", 0, 1.0, "oxygen"]
  ]}
```
The response is as follows...
```{ <reaction_id> : [ [reaction], False, True, "brought into balance"  ]  }```
...where the list elements contain the reaction, whether or not it was originally balanced, whether or not it is now balanced, and a message from the balancer.

For example...
```
{ "photosynthesis" :
  [
    [	["H2O", 0, -6.0, "water"],
      ["CO2", 0, -6.0, "carbon dioxide"],
      ["C6H12O6", 0, 1.0, "glucose"],
      ["O2", 0, 6.0, "oxygen"]
    ],
   False,
   True,
  "brought into balance"
  ]
}

```

### SBML
The SBML endpoint at ```http://<URL>:<PORT>/balance/SBML``` expects SBML Level 3 with the 'fbc' package used to store molecular species' ```chemicalFormula``` and ```charge```. Earlier SBML levels will be rejected by the service. Use POST to send SBML as data.

Existing balanced reactions maintain original stoichiometries and tagged (in ```Notes```) as 'already balanced'.
Unbalanced reactions brought into balance are tagged as 'brought into balance'.
Unbalanced reactions not brought into balance are tagged as 'not brought into balance'.
Reactions that could not be processed because one or more components lacked formula or charge are tagged 'could not be assessed'.
