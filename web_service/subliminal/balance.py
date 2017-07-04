'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import itertools, chem_utils, math_utils


def balance_model(model, verbose=True):
    '''Attempts to mass / charge balance models.'''
    for reaction in model.reactions:
        # If NOT exchange reaction:
        if len(reaction.metabolites) > 1 and \
                len(check_reaction_balance(reaction)) > 1:
            reac_def = [(met.formula, met.charge, stoich, met.id)
                        for met, stoich in reaction.metabolites.iteritems()]
            result = balance_reac(reac_def, max_stoich=10)

            # If reaction has been 'fixed' update stoichiometries:
            if result[0] and not result[1]:
                if verbose:
                    print 'Reaction %s fixed:' % reaction.id
                    print 'FROM: ' + reaction.build_reaction_string()

                metabolites = {met: stoich
                               for met, stoich
                               in reaction.metabolites.iteritems()}
                reaction.clear_metabolites()

                for val in result[2]:
                    metabolite = [met for met in metabolites
                                  if met.id == val[3]][0]
                    metabolites[metabolite] = val[2]

                reaction.add_metabolites(metabolites)

                if verbose:
                    print 'TO:   ' + reaction.build_reaction_string()
            elif not result[0]:
                if verbose:
                    print 'Reaction %s unbalanced' % reaction.id

    return model


def balance_reac(reaction_def, optional_comp=None, max_stoich=8.0):
    '''Applies linear programming to balance reaction.'''
    if optional_comp is None:
        optional_comp = [('H', 1, 'proton_BALANCER'), ('H2O', 0, 'water_BALANCER')]

    all_formulae = [[x[0] for x in reaction_def if x[2] <= 0],
                    [x[0] for x in reaction_def if x[2] > 0]]
    all_formulae.extend([[x[0] for x in optional_comp] for _ in range(2)])

    all_charges = [[x[1] for x in reaction_def if x[2] <= 0],
                   [x[1] for x in reaction_def if x[2] > 0]]
    all_charges.extend([[x[1] for x in optional_comp] for _ in range(2)])

    all_ids = [[x[3] for x in reaction_def if x[2] <= 0],
               [x[3] for x in reaction_def if x[2] > 0]]
    all_ids.extend([[x[2] for x in optional_comp] for _ in range(2)])

    all_elem_comp = [[_get_elem_comp(formula, idx)
                      for formula in formulae]
                     for idx, formulae in enumerate(all_formulae)]

    (balanced, stoichs) = _optimise(_get_elem_matrix(all_elem_comp,
                                                     all_charges),
                                    max_stoich, len(optional_comp))

    balanced_def = sorted(_get_reaction_def(stoichs, all_formulae,
                                            all_charges, all_ids)) \
        if balanced else reaction_def

    return balanced, \
        balanced and _compare_reaction_defs(reaction_def, balanced_def), \
        balanced_def


def check_model_balance(model):
    '''Checks a model for mass / charge balancing.'''
    model_balance = {}

    for reaction in model.reactions:
        # If NOT exchange reaction:
        if len(reaction.metabolites) > 1:
            # Check mass balance and remove tiny coefficients due to rounding:
            reaction_balance = check_reaction_balance(reaction)

            if len(reaction_balance) > 0:
                model_balance[reaction.id] = reaction_balance

    return model_balance


def check_reaction_balance(reaction):
    '''Checks a reaction for mass / charge balancing.'''
    if len(reaction.metabolites) > 1:
        # Check mass balance and remove tiny coefficients due to rounding:
        mass_balance = reaction.check_mass_balance()

        return {key: val for key, val in mass_balance.iteritems()
                if abs(val) > 1e-6}


def _get_reaction_def(stoichs, all_formulae, all_charges, all_ids):
    '''Formats the input into (formula, charge, stoichiometry) tuples.'''
    return [(a[0], b, a[1] * c, d)
            for a, b, c, d in zip([(x, -1 if idx % 2 == 0 else 1)
                                   for idx, formulae in enumerate(all_formulae)
                                   for x in formulae],
                                  [y for charges in all_charges
                                   for y in charges],
                                  _simplify_stoichs(stoichs),
                                  [z for ids in all_ids for z in ids])
            if c > 1e-8]


def _compare_reaction_defs(def1, def2):
    '''Returns True/False depending on whether reaction definitions are
    equal.'''
    return def1 is not None and def2 is not None and \
        len(def1) == len(def2) and \
        all([x[0] == y[0] and x[1] == y[1] and math_utils.isclose(x[2], y[2])
             for x, y in zip(sorted(def1), sorted(def2))])


def _simplify_stoichs(stoichs):
    '''Attempts to simplify stoichs of 1.00000001 to 1.0.'''
    return [float(round(stoich)) if math_utils.isclose(stoich, round(stoich))
            else stoich for stoich in stoichs]


def _get_elem_comp(formula, idx):
    '''Returns elemental compositions, multiplying stoichiometry by -1
    if reactants.'''
    elem_comp = chem_utils.get_elem_comp(formula)
    elem_comp.update((x, y * (-1 if idx % 2 == 0 else 1))
                     for x, y in elem_comp.items())
    return elem_comp


def _get_elem_matrix(all_elem_comp, all_charges):
    '''Gets the elemental (and charge) matrix representing rows of elements
    and columns of compounds.'''
    a_matrix = []
    elements = [elem_comp.keys() for elem_comps in all_elem_comp
                for elem_comp in elem_comps]

    for element in set([item for sublist in elements for item in sublist]):
        a_matrix.append([elem_comp[element] if element in elem_comp else 0
                         for elem_comps in all_elem_comp
                         for elem_comp in elem_comps])

    a_matrix.append(list(itertools.chain(*([charge * (-1 if idx % 2 == 0
                                                      else 1)
                                            for charge in charges]
                                           for idx, charges
                                           in enumerate(all_charges)))))

    return a_matrix


def _optimise(a_matrix, max_stoich, num_opt_formulae):
    '''Optimise linear program and return result.'''
    bounds = [(1.0, max_stoich)] * (len(a_matrix[0]) -
                                    num_opt_formulae * 2) + \
        [(0.0, max_stoich)] * num_opt_formulae * 2

    return math_utils.linprog([1] * len(a_matrix[0]),
                              a_matrix,
                              [0.0] * (len(a_matrix)),
                              bounds)
