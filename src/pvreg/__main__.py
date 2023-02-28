'''
pvreg
ver 1.2.0
2023.02.28
everythingthatcounts@gmail.com
'''

import numpy as np
import pandas as pd
import multiprocessing as mp
import statsmodels.api as sm
import enlighten
import re
import copy
import json
import argparse

# Classes

class Group:
    '''Group information'''
    def __init__(self, name, value, mean, sd):
        self.name = name # name of the grouping variable
        self.value = value # value of the grouping variable inticating class instance
        self.mean = mean # mean of a priori distribution of theta
        self.sd = sd # sd of ^^^^

class Item:
    '''Item information'''
    def __init__(self, name, pars, par_labs, model, resp_cats=None):
        self.name = name # name of the item variable
        self.model = model # 2plm,3plm,pcm,gpcm
        self.pars = pars # current estimates of item parameters [a,b]|[a,b,c]|[a,b1,...,bmax]
        self.par_labs = par_labs # string labels for item parameters
        if resp_cats:
            self.resp_cats = resp_cats # values of item responses (default: [0,1,..,max])
            self.n_cats = len(resp_cats) # number of response categories
        else:
            if model in ['1plm','2plm', '3plm']:
                self.n_cats=2
            else:
                self.n_cats=len(pars)
            self.resp_cats = list(range(self.n_cats))

# Functions

def make_quad(group, nip):
    '''GH normal quadrature points and probabilities for specified number of integration points '''
    GH_quad_pts, GH_quad_probs = np.polynomial.hermite.hermgauss(nip)
    x = GH_quad_pts * np.sqrt(2) * group.sd + group.mean # quad points
    p = GH_quad_probs / np.sqrt(np.pi) # quad probabilities
    return x, p


def item_prob(item, resp, theta):
    '''Probability of responding resp to an item when ability is theta'''
    if item.n_cats == 2:
        if item.model=='3plm':
            p = item.pars[2] + (1. - item.pars[2]) / (1. + np.exp(-item.pars[0] * (theta - item.pars[1])))
            if resp==item.resp_cats[1]:
                return p
            else:
                return 1. - p
        else:
            p = 1. / (1. + np.exp(-item.pars[0] * (theta - item.pars[1])))
            if resp==item.resp_cats[1]:
                return p
            else:
                return 1. - p
    else:
        if item.model == 'grm':
            if resp==item.resp_cats[0]:
                return 1. - 1. / (1. + np.exp(-item.pars[0] * (theta - item.pars[1])))
            elif  resp==item.resp_cats[-1]:
                return  1. / (1. + np.exp(-item.pars[0] * (theta - item.pars[-1])))
            else:
                pos = item.resp_cats.index(resp)
                return 1. / (1. + np.exp(-item.pars[0] * (theta - item.pars[pos]))) - 1. / (1. + np.exp(-item.pars[0] * (theta - item.pars[pos+1])))
        else: # gpcm
            expsum_all = 1.
            for pos in range(1,item.n_cats):
                expsum_all += np.exp(item.pars[0] * ( pos*theta - sum(item.pars[1:pos+1]) ))
            if resp==item.resp_cats[0]:
                return  1. / expsum_all
            else:
                pos = item.resp_cats.index(resp)
                return np.exp(item.pars[0] * ( pos*theta - sum(item.pars[1:pos+1]) )) / expsum_all


def load_data(in_estimates_df, in_responses_df,
              out_estimates_df, out_responses_df,
              student_exog_df,
              student_id, school_id, fixed_effects):
    ''' Returns:
        (1) dict with Group instances (key=group.value)
        (2) dict with Item instances (key=item.name)
        (3) dict from item responses and grouping variable; key=(item.name,item_response,group.value), value = index
        (4) dataframe with combined IRT-parameters estimates and covariance matrix
        (5) dataframe with [student_id, school_id, 'group','exam_var']+fixed_efects+fixed_effects_X_exam columns
        (6) the fixed_efects+fixed_effects_X_exam list'''

    # combining estimates and responses (recoding group values if repeated)
    g_pars_in = in_estimates_df.loc[in_estimates_df['par'].str.contains('_theta'), ['var', 'par', 'est']]
    group_var_in = re.search("[^\b]*(?=_\d*$)", g_pars_in.iloc[0, 0])[0]
    g_values_in = set(g_pars_in['var'].apply(lambda x: int(re.search("(\d*$)", x)[0])))
    g_pars_out = out_estimates_df.loc[out_estimates_df['par'].str.contains('_theta'), ['var', 'par', 'est']]
    group_var_out = re.search("[^\b]*(?=_\d*$)", g_pars_out.iloc[0, 0])[0]
    g_values_out = set(g_pars_out['var'].apply(lambda x: int(re.search("(\d*$)", x)[0])))
    common_g_values = g_values_in.intersection(g_values_out)
    if len(common_g_values):
        new_val = max(g_values_in.union(g_values_out)) + 1
        g_out_recode = {}
        for old_val in common_g_values:
            g_out_recode[old_val] = new_val
            new_val += 1
        for old_val, new_val in g_out_recode.items():
            out_estimates_df.loc[out_estimates_df.loc[:, 'var'] == group_var_out+"_" + str(old_val), 'var'] = "group_" + str(
                new_val)
            out_estimates_df.rename(
                columns={group_var_out+"_" + str(old_val) + "_mean_theta": "group_" + str(new_val) + "_mean_theta",
                         group_var_out+"_" + str(old_val) + "_sd_theta": "group_" + str(new_val) + "_sd_theta"}, inplace=True)
            out_responses_df.loc[out_responses_df[group_var_out]==old_val, group_var_out] = new_val

    #unifying group name for in_estimates_df (already done for out_estimates_df)
    for g in g_values_in:
        in_estimates_df.loc[
            in_estimates_df.loc[:, 'var'] == group_var_in + "_" + str(g), 'var'] = "group_" + str(
            g)
        in_estimates_df.rename(
            columns={group_var_in + "_" + str(g) + "_mean_theta": "group_" + str(g) + "_mean_theta",
                     group_var_in + "_" + str(g) + "_sd_theta": "group_" + str(g) + "_sd_theta"},
            inplace=True)

    estimates_df = pd.concat([in_estimates_df, out_estimates_df])
    cols_tofillna = [i for i in estimates_df.columns if i not in ['var', 'par', 'est', 'groupN_itemCATS']]
    estimates_df[cols_tofillna] = estimates_df[cols_tofillna].fillna(0.)

    #unifying group name for responses
    in_responses_df.rename(
                columns={group_var_in: "group"}, inplace=True)
    out_responses_df.rename(
                columns={group_var_out: "group"}, inplace=True)

    in_responses_df['exam_var'] = 0  # I wish so much no item is called 'exam_var', also - changes external df
    out_responses_df['exam_var'] = 1  # I wish so much no item is called 'exam_var', also - changes external df
    responses_df = pd.concat([in_responses_df, out_responses_df]).reset_index(drop=True)

    # group instances
    groups = {}
    g_pars = estimates_df.loc[estimates_df['par'].str.contains('_theta'), ['var', 'par', 'est']]
    g_values = set(g_pars['var'].apply(lambda x: int(re.search("(\d*$)", x)[0])))
    for g in g_values:
        g_pars_g = g_pars.loc[(g_pars['var'] == 'group_' + str(g)), ['par', 'est']].set_index('par')
        groups[g] = Group("group", g, g_pars_g.loc['mean_theta'].values[0], g_pars_g.loc['sd_theta'].values[0])

    # item instances
    items = {}
    i_pars = estimates_df.loc[~estimates_df['par'].str.contains('_theta'), ['var', 'par', 'est','groupN_itemCATS']]
    i_names = set(i_pars['var'])
    for i in i_names:
        i_pars_i = i_pars.loc[(i_pars['var'] == i), ['par', 'est']]
        parameters_i = list(i_pars_i['est'])
        parlabels_i = list(i_pars_i['par'].apply(lambda x: re.search("(?<=_)\w*", x)[0]))
        model_i = list(i_pars_i['par'].apply(lambda x: re.search("\w*(?=_)", x)[0]))[0]
        cats_i = [int(i) for i in list(i_pars.loc[i_pars['var']==i, 'groupN_itemCATS'].dropna())]
        items[i] = Item(i, parameters_i, parlabels_i, model_i, cats_i)

    # plugging in the item response indices into groups
    responses = {}
    for g in groups.values():
        g.indx = responses_df.loc[responses_df["group"] == g.value, :].index.tolist()
        g.n = len(g.indx)  # total number of observations for the group
        g_items = {}
        for item_name, item_data in responses_df.iloc[g.indx,
                                                      ~responses_df.columns.isin(
                                                          [student_id, 'group', 'exam_var'])].iteritems():
            i_n = item_data.count()
            if i_n:
                g_items[item_name] = i_n
                i = items[item_name]
                for r in i.resp_cats:
                    responses[(i.name, r, g.value)] = item_data.loc[item_data == r].index.tolist()
        g.items = g_items

    # collecting student data
    student_df = responses_df.loc[:, [student_id, 'group', 'exam_var']]
    student_df = student_df.merge(student_exog_df[[student_id,school_id]+fixed_effects], how='left', on=student_id).sort_index()
    # dealing with missing data
    for fixed_effect in fixed_effects:
        if student_df[fixed_effect].isnull().sum():
            fe_mean=student_df[fixed_effect].mean()
            student_df[fixed_effect] = student_df[fixed_effect].fillna(
                student_df.groupby(school_id)[fixed_effect].transform('mean')) # missing values are filled with school means
            if student_df[fixed_effect].isnull().sum():
                student_df[fixed_effect] = student_df[fixed_effect].fillna(fe_mean)  # if that does not help - global mean
    # interaction of exog with measurement occasion
    for fixed_effect in fixed_effects:
        student_df[fixed_effect+'_X_exam']=student_df[fixed_effect]*student_df['exam_var']
    fixed_effects_all = fixed_effects + [fixed_effect+'_X_exam' for fixed_effect in fixed_effects]

    return groups, items, responses, estimates_df, student_df, fixed_effects_all


def likelihood_long(items, groups, responses, theta):
    '''Likelihood of a response vectors for given thetas'''

    l = np.ones_like(theta, dtype=float)

    for g in groups.values():
        for item in g.items.keys():
            i = items[item]
            for r in i.resp_cats:
                indx = responses[(i.name, r, g.value)]
                if len(indx):
                    l[indx] *= item_prob(i, r, theta[indx])
    return l


def likelihood_wide(items, groups, responses, nip):
    '''Likelihood of response vectors at given quadrature_points'''

    l = np.ones((sum([g.n for g in groups.values()]), nip))

    for g in groups.values():
        x = make_quad(g, nip)[0]
        for item in g.items.keys():
            i = items[item]
            for r in i.resp_cats:
                indx = responses[(i.name, r, g.value)]
                if len(indx):
                    l[indx, :] *= item_prob(i, r, x)
    return l


def eap(items, groups, responses, nip):
    '''Expected a posteriori estimate and its standard error'''

    l = likelihood_wide(items, groups, responses, nip)

    n = sum([g.n for g in groups.values()])
    est_theta = np.empty(n)
    est_theta_se = np.empty(n)

    for g in groups.values():
        x, p = make_quad(g, nip)

        prob_g = l[g.indx, :]
        prob_g *= p
        prob_g /= np.reshape(np.repeat(np.sum(prob_g, axis=1), nip), (g.n, nip))  # TDL: maybe smthng more elegant?

        est_theta[g.indx] = np.inner(prob_g, x)

        est_theta_se[g.indx] = np.sqrt(np.inner(prob_g, x ** 2) - est_theta[g.indx] ** 2)

    return est_theta, est_theta_se


def normal_den(x, mu, sigma): #20x faster than: from scipy.stats import norm
    '''Probability of N(x,mu,sigma), vector imputs handled rowwise'''
    return np.exp( - ( np.square((x-mu)) / ( 2*np.square(sigma) ) ) ) / ( sigma * np.sqrt(2*np.pi) )


def draw_mcmc(items, groups, responses, theta_t, proposal_sd, prior_mean, prior_sd):
    theta_tt = np.random.normal(theta_t, proposal_sd)

    L_t = likelihood_long(items, groups, responses, theta_t) * normal_den(theta_t, prior_mean, prior_sd)
    L_tt = likelihood_long(items, groups, responses, theta_tt) * normal_den(theta_tt, prior_mean, prior_sd)

    alpha = np.minimum(1, L_tt / L_t)

    updates = np.arange(len(theta_tt))[np.random.random(len(theta_tt)) > alpha]

    if len(updates):
        theta_tt[updates] = theta_t[updates]

    return theta_tt


def noisify_estimates(items, groups, estimates):
    '''Replaces estimates of groups and items with random draw drop N(est,cov_est)'''

    # absolutely no bother, that this 'noisy_est' column is not dropped afterwards (reused and overwritten many times)
    estimates['noisy_est'] = np.random.multivariate_normal(estimates['est'],
                                                           estimates.loc[:,
                                                           list(estimates['var'] + '_' + estimates['par'])])

    g_pars = estimates.loc[estimates['par'].str.contains('_theta'), ['var', 'par', 'noisy_est']]
    for g in groups.keys():
        g_pars_g = g_pars.loc[(g_pars['var'] == 'group_' + str(g)), ['par', 'noisy_est']].set_index('par')

        groups[g].mean = g_pars_g.loc['mean_theta'].values[0]

        sd_g = g_pars_g.loc['sd_theta'].values[0]
        # out-of-parameter-space-boundary corrections:
        if sd_g <= 0:  # we do not want negative sd (can't imagine a scenario where we would have one fixed at zero as well)
            sd_g = 0.0001

        groups[g].sd = sd_g

    i_pars = estimates.loc[~estimates['par'].str.contains('_theta'), ['var', 'par', 'noisy_est']]
    for i in items.keys():
        i_pars_i = i_pars.loc[(i_pars['var'] == i), ['par', 'noisy_est']]
        parameters_i = list(i_pars_i['noisy_est'])

        # out-of-parameter-space-boundary corrections:
        if items[
            i].model == '3plm':  # for 3plm we want the c parameter to be within (0,1) range; allowing 0 just in case of a lazy 2plm
            if parameters_i[2] < 0:
                parameters_i[2] = 0.0001
            if parameters_i[2] >= 1:
                parameters_i[2] = 0.9999
        # Comment: (1) Correcting instances of a<0 is skipped to allow for a lazy multigroup with selection variables as items
        #         (2) Correcting for misordering of difficulties in grm is on TLD (since forever)

        items[i].pars = parameters_i

def reff_from_multilevel_results(model,results,school_id, exam_var):
    reff = results.random_effects
    reff_c = results.random_effects_cov
    reff_school = {}
    reff_school_cov = {}
    for key in model.group_labels:
        reff_school[key] = reff[key][[school_id, exam_var]]
        reff_school_cov[key] = reff_c[key].loc[[school_id, exam_var], [school_id, exam_var]]
    return reff_school, reff_school_cov

def priors_from_multilevel_results(model, results, data, school_id, student_id, exam_var):
    ''' (1) Converts MixedLM estimates into a useful Pandas Dataframe
        (2) Adds prior_mean and its randomized version (a draw according to standard errors of effect estimates);
        (3) Adds prior_sd - sqrt of sum of error variances'''

    reff = results.random_effects
    reff_c = results.random_effects_cov

    # random effect predictions
    n_reff = sum([len(reff[key]) - 2 for key in model.group_labels])
    df_reff = np.zeros((n_reff, 11))

    r_0 = 0
    r_n = 0
    for key in model.group_labels:

        r_n += len(reff[key]) - 2
        slice_key = range(r_0, r_n)

        student_labels = reff[key].index[reff[key].index.str.contains('\[')].tolist()

        # school_id and student_id values
        df_reff[slice_key, 0] = key
        df_reff[slice_key, 1] = [int(re.search('(?<=\[)\d*(?=\])', label)[0]) for label in student_labels]

        # random effects estimates
        df_reff[slice_key, 2] = reff[key][school_id]
        df_reff[slice_key, 3] = reff[key][exam_var]
        df_reff[slice_key, 4] = reff[key][student_labels]

        # variances of random effects estimates (turned to SE, at the end)
        df_reff[slice_key, 5] = reff_c[key].loc[school_id, school_id]
        df_reff[slice_key, 6] = reff_c[key].loc[exam_var, exam_var]
        df_reff[slice_key, 7] = np.diag(reff_c[key].loc[student_labels, student_labels])

        # randomized random effects estimates
        try:
            reff_rand = pd.Series(np.random.multivariate_normal(reff[key], reff_c[key], check_valid='ignore'), index=reff[key].index)
        except:
            stable_cov = copy.deepcopy(reff_c[key])
            stable_cov.loc[student_labels, student_labels] = np.array(
                stable_cov.loc[student_labels, student_labels]).diagonal() * np.eye(len(student_labels))
            stable_cov.loc[[school_id,exam_var], student_labels] = 0
            stable_cov.loc[student_labels,[school_id,exam_var]] = 0
            stable_cov=(stable_cov+stable_cov.T)/2
            reff_rand = pd.Series(np.random.multivariate_normal(reff[key], stable_cov), index=reff[key].index)
        df_reff[slice_key, 8] = reff_rand[school_id]
        df_reff[slice_key, 9] = reff_rand[exam_var]
        df_reff[slice_key, 10] = reff_rand[student_labels]

        r_0 = r_n

    colnames = [school_id, student_id,
                school_id + '_reff', exam_var + '_reff', student_id + '_reff',
                school_id + '_reff_se', exam_var + '_reff_se', student_id + '_reff_se',
                school_id + '_reff_rand', exam_var + '_reff_rand', student_id + '_reff_rand',
                ]

    df_reff = pd.DataFrame(df_reff, columns=colnames)

    # fixed effects predictions
    b = np.array(results.fe_params)
    b_names = results.fe_params.index.tolist()
    b_cov = results.cov_params().loc[b_names, b_names]

    x = np.concatenate((np.ones((len(data), 1)),
                        data.loc[:, b_names[1:]].to_numpy()),
                       axis=1)  # adding intercept (always first)
    x_uniq = np.unique(x, axis=0)

    xb_uniq = x_uniq.dot(b)
    xb_ve_uniq = [obs.dot(b_cov.dot(obs)) for obs in x_uniq]
    xb_rand_uniq = np.random.normal(xb_uniq, np.sqrt(xb_ve_uniq))

    df_fixed = np.zeros((len(data), 3))
    for i, u in enumerate(x_uniq):
        indx = np.where((x == u).all(axis=1))
        df_fixed[indx, :] = [xb_uniq[i], xb_ve_uniq[i], xb_rand_uniq[i]]

    # merging with original data index
    df = data.loc[:, [school_id, student_id]]
    df[['xb', 'xb_se', 'xb_rand']] = df_fixed
    df = df.merge(df_reff, how='left', on=[school_id, student_id])

    # exam_var==0 cases have 0 prediction on blup (error is dropped as well)
    zeros_index = list(data[exam_var] == 0)
    df.loc[zeros_index, [exam_var + '_reff', exam_var + '_reff_se', exam_var + '_reff_rand']] = np.zeros(
        (sum(zeros_index), 3))

    # adding priors
    df['prior_mean'] = df[['xb', school_id + '_reff', exam_var + '_reff', student_id + '_reff']].sum(axis=1)
    df['prior_mean_rand'] = df[
        ['xb_rand', school_id + '_reff_rand', exam_var + '_reff_rand', student_id + '_reff_rand']].sum(axis=1)
    df['prior_sd'] = np.sqrt(
        df[['xb_se', school_id + '_reff_se', exam_var + '_reff_se', student_id + '_reff_se']]
        .sum(axis=1))
    df[['xb_se', school_id + '_reff_se', exam_var + '_reff_se', student_id + '_reff_se']] = np.sqrt(df[[
        'xb_se', school_id + '_reff_se', exam_var + '_reff_se', student_id + '_reff_se']])

    #instead of the whole df I return only priors - rest is unused
    return df['prior_mean_rand'], df['prior_sd']

def multilevel(data, theta_var, school_id, student_id, exam_var, fixed_effects=[]):
    '''fits MixedLM according to 3-level school EVA model'''
    fixed = '1 + ' + exam_var
    for eff in fixed_effects:
        fixed += ' + ' + eff

    vc = {student_id: '0+C(' + student_id + ')'}

    model = sm.MixedLM.from_formula(theta_var + ' ~ ' + fixed,
                                    groups=school_id,
                                    re_formula='1 + ' + exam_var,
                                    vc_formula=vc,
                                    data=data)
    results = model.fit()

    return model, results


def rescale_priors(prior_mean, prior_sd, groups):
    for g in groups.values():
        _mean = prior_mean[g.indx].mean()
        _v = (prior_sd[g.indx] ** 2 + prior_mean[g.indx] ** 2).mean() - _mean ** 2
        _sd = np.sqrt(_v)
        prior_mean[g.indx] = (prior_mean[g.indx] - _mean) * (g.sd / _sd) + g.mean
        prior_sd[g.indx] = prior_sd[g.indx] * (g.sd / _sd)

    return prior_mean, prior_sd


def child_initialize(_items, _groups, _responses, _estimates, _burn0, _burn1,
                        _draw_delta,
                        _keep_pv, _keep_reff, _mlv_data, _mlv_args):
    '''initializer for mp.Pool()'''
    global items, groups, responses, estimates, burn0, burn1, draw_delta, keep_pv, keep_reff, mlv_data, mlv_args
    items = _items
    groups = _groups
    responses= _responses
    estimates = _estimates
    burn0 = _burn0
    burn1 = _burn1
    draw_delta = _draw_delta
    keep_pv = _keep_pv
    keep_reff = _keep_reff
    mlv_data = _mlv_data
    mlv_args = _mlv_args

def generate_pv_child(n_draw, max_indep_chains, shift):
    '''2-argument version of generate_pv_single() for pool.apply()'''
    return generate_pv_single(items, groups, responses, estimates, burn0, burn1,
                                            n_draw, draw_delta, max_indep_chains,
                                            keep_pv, keep_reff, mlv_data, mlv_args,
                                            shift)
def generate_pv_multipro(items, groups, responses, estimates, njobs=2,
                        burn0=10,burn1=20, n_draw=10, draw_delta=10, max_indep_chains=10,
                        keep_pv=True,
                        keep_reff=True, mlv_data=[], mlv_args=[]):
    '''Generates MCMC chains for regression conditioned PVs in parallel mode (njobs)'''

    if max_indep_chains < njobs:
        max_indep_chains = njobs

    if n_draw < njobs:
        njobs = n_draw

    # list: [n_draw, max_indep_chains, shift]
    pool_args = [[(job + n_draw) // njobs, (job + max_indep_chains) // njobs] for job in range(njobs)]
    pool_args.reverse()
    pool_args = [pool_args[job]
                         + [sum([pool_args[j][0] for j in range(job)])]
                         for job in range(njobs)]


    pool = mp.Pool(njobs, initializer = child_initialize,
                   initargs = (items, groups, responses, estimates, burn0, burn1,
                               draw_delta, keep_pv, keep_reff, mlv_data, mlv_args)
                   )
    pool_res =   [pool.apply_async(generate_pv_child, args=(n_draw, max_indep_chains, shift))
                    for n_draw, max_indep_chains, shift in pool_args]
    pool_res = [x.get() for x in pool_res]
    pool.close()

    pvs = {}
    reffs_school = {}
    reffs_school_cov = {}
    for i in range(len(pool_res)):
        pvs.update( pool_res[i][0] )
        reffs_school.update( pool_res[i][1] )
        reffs_school_cov.update( pool_res[i][2] )

    return pvs, reffs_school, reffs_school_cov


def generate_pv_single(items, groups, responses, estimates,
                burn0=10, burn1=20, n_draw=10, draw_delta=10, max_indep_chains=10,
                keep_pv=True,
                keep_reff=True, mlv_data=[], mlv_args=[],
                shift=-1):
    '''Generates MCMC chains for regression conditioned PVs '''

    pvs={}
    reffs_school = {}
    reffs_school_cov = {}

    if n_draw > max_indep_chains:
        max_chain = max_indep_chains
    else:
        max_chain = n_draw
    chain_draws = [list(range(draw_delta, 1+draw_delta * (n_draw // max_chain + (n_draw % max_chain > chain)), draw_delta))
               for chain in range(max_chain)]

    #progress barr
    if shift in [-1,0]:
        manager = enlighten.get_manager()
        bar_desc="MCMC"+(shift+1)*" (job 1)"
        progress_bar = manager.counter(total=sum([burn0+burn1 + max(i) + 1 for i in chain_draws]),
                                       desc=bar_desc,
                                       unit="iter")
        shift=0

    for chain in range(max_chain):

        # randomising parameter estimates according to their covariance matrix
        _groups = copy.deepcopy(groups)
        _items = copy.deepcopy(items)
        noisify_estimates(_items, _groups, estimates)

        theta_t, proposal_sd = eap(_items, _groups, responses, 91)  # paranoid nip - TDL scale it down at some time

        prior_mean = np.empty_like(theta_t)
        prior_sd = np.empty_like(theta_t)
        for g in _groups.values():
            prior_mean[g.indx] = g.mean
            prior_sd[g.indx] = g.sd

        # burn
        for _ in range(burn0):
            theta_t = draw_mcmc(_items, _groups, responses, theta_t, proposal_sd, prior_mean, prior_sd)
            if shift==0:
                progress_bar.update()

        # multilevel burn
        if len(mlv_args):
            school_id, student_id, exam_var, fixed_effects = mlv_args
            for _ in range(burn1):
                mlv_data['theta'] =theta_t
                mlv_model, mlv_results = multilevel(mlv_data, 'theta', *mlv_args)
                prior_mean, prior_sd = priors_from_multilevel_results(mlv_model, mlv_results, mlv_data,
                                                                           school_id, student_id, exam_var)
                prior_mean, prior_sd = rescale_priors(prior_mean, prior_sd,
                                                         _groups)  #TLD - make sure does its job
                theta_t = draw_mcmc(_items, _groups, responses, theta_t, proposal_sd, prior_mean, prior_sd)
                if shift == 0:
                    progress_bar.update()

        # final chain
        for draw in range(max(chain_draws[chain]) + 1):

            theta_t = draw_mcmc(_items, _groups, responses, theta_t, proposal_sd, prior_mean, prior_sd)

            if draw in chain_draws[chain]:

                if len(mlv_args):
                    mlv_data['theta'] = theta_t
                    mlv_model, mlv_results = multilevel(mlv_data, 'theta', *mlv_args)
                    prior_mean, prior_sd = priors_from_multilevel_results(mlv_model, mlv_results, mlv_data,
                                                                          school_id, student_id, exam_var)
                    prior_mean, prior_sd = rescale_priors(prior_mean, prior_sd,
                                                          _groups)

                col = shift + sum([len(chain_draws[d]) for d in range(chain)]) + chain_draws[chain].index(draw)

                if keep_pv:
                    pvs[col] = theta_t

                if keep_reff:
                    reffs_school[col], reffs_school_cov[col] = reff_from_multilevel_results(mlv_model, mlv_results,school_id, exam_var)
            if shift == 0:
                progress_bar.update()

    return pvs, reffs_school, reffs_school_cov


def pvs_to_reff_df(pv_results, school_id):
    '''Transforms pv_results into a dataframe with pvs and a dataframe with school effect estimates'''
    pvs = pv_results[0]
    reffs_school = pv_results[1]
    reffs_school_cov = pv_results[2]
    m = len(reffs_school)
    reff_df = []
    for schl in reffs_school[0].keys():
        b = np.array([reffs_school[pv][schl] for pv in range(m)])
        b_est = np.mean(b, axis=0)
        W = np.zeros(shape=(len(b_est), len(b_est)))
        B = np.zeros_like(W)
        for pv in range(m):
            W += reffs_school_cov[pv][schl]
            B += np.outer(b[pv] - b_est, b[pv] - b_est)
        W /= m
        B /= m - 1
        V = W + ((1 + m) / m) * B
        b_est_se = np.sqrt(np.diag(V))
        reff_df.append([schl,
                        b_est[0],
                        b_est_se[0],
                        b_est[1],
                        b_est_se[1],
                        V.values.tolist(),
                        W.values.tolist(),
                        B.tolist(),
                        m
                        ])
    reff_df = pd.DataFrame(reff_df,
                          columns=[school_id, 'mean_schl', 'mean_schl_se', 'eva_schl', 'eva_schl_se', 'COV', 'W', 'B',
                                    'm'])
    pvs=pd.DataFrame(pvs)

    return reff_df, pvs

def pvreg(estimates_in, responses_in, # dataframes with uirt estimates and item responses - exam 0
          estimates_out, responses_out, # dataframes with uirt estimates and item responses - exam 1
          student_exog, # dataframe with student variables and school_id
          student_id, # str - name of student id column (common in responses_# and student_exog df)
          school_id, # str - name of school id column (in student_exog df)
          fixed_effects=[], # list with column names of exog variables to be included in fixed effects part (in student_exog df)
          npv = 10, # number of plausible values to be drawn
          keep_pv =True, # whether to keep plausible values
          csvname = '', # prefix used for saved csv files (can include a path)
          njobs = 1
         ):

    groups, items, responses, estimates_df, student_df, fixed_effects_all = load_data(estimates_in, responses_in,
                                                                   estimates_out, responses_out,
                                                                   student_exog,
                                                                   student_id, school_id, fixed_effects)
    if njobs==1:
        pv_results = generate_pv_single(items, groups, responses, estimates_df,
                                    burn0=10, burn1=20, n_draw=npv, draw_delta=10, max_indep_chains=npv,
                                    keep_pv=keep_pv,
                                    keep_reff=True, mlv_data=student_df,
                                    mlv_args=[school_id, student_id, 'exam_var', fixed_effects_all],
                                    shift=-1)
    else:
        pv_results = generate_pv_multipro(items, groups, responses, estimates_df, njobs=njobs,
                             burn0=10, burn1=20, n_draw=npv, draw_delta=10, max_indep_chains=npv,
                             keep_pv=keep_pv,
                             keep_reff=True, mlv_data=student_df,
                             mlv_args=[school_id, student_id, 'exam_var', fixed_effects_all])

    reffs, pvs = pvs_to_reff_df(pv_results, school_id)

    reffs.to_csv(csvname+'reffs.csv', index=False,
                 columns=[school_id, 'mean_schl', 'mean_schl_se', 'eva_schl', 'eva_schl_se'
                     ,'COV', 'W', 'B', 'm'])
    if keep_pv:
        pvs.columns = ['pv_'+str(i) for i in pvs.columns]
        pv_cols=[student_id, school_id, 'exam_var'] + pvs.columns.to_list()
        pvs[student_id] = student_df[student_id]
        pvs[school_id] = student_df[school_id] #unnecesary but useful
        pvs['exam_var'] = student_df['exam_var']
        pvs[pv_cols].to_csv(csvname+'pvs.csv', index=False)

 #   return reffs

def run_pv_reg(args):

    file_path = args.path

    estimates_in_df = pd.read_csv(file_path+args.estimates_in)
    responses_in_df = pd.read_csv(file_path+args.responses_in)

    estimates_out_df = pd.read_csv(file_path+args.estimates_out)
    responses_out_df = pd.read_csv(file_path+args.responses_out)

    student_exog_df = pd.read_csv(file_path+args.student_exog)

    pvreg(
        estimates_in = estimates_in_df,
        responses_in = responses_in_df,
        estimates_out = estimates_out_df,
        responses_out = responses_out_df,
        student_exog = student_exog_df,
        student_id = args.student_id,
        school_id = args.school_id,
        fixed_effects = args.fixed_effects.split(),
        npv = int(args.npv),
        keep_pv = args.keep_pv == 'True',
        csvname = args.out_path+args.out_files,
        njobs = int(args.njobs)
    )

def main():
    parser = argparse.ArgumentParser(prog='pvreg',
                                     description='Generate PVs conditioned on multilevel regression')

    parser.add_argument('-c', '--conf_file',
                        type=str,
                        help='JSON config file',
                        metavar="FILE")

    parser.add_argument('-p', '--path',
                        type=str,
                        help='Path where input data is stored',
                        required=False)
    parser.add_argument('--out_path',
                        type=str,
                        help='Path where output is stored (default = args.path)',
                        required=False)
    parser.add_argument('--njobs',
                        default=1, type=int,
                        help='Number of parallel jobs for multiprocessing',
                        required=False)
    parser.add_argument('--npv',
                        default=15, type=int,
                        help='Number of plausible values drawn',
                        required=False)
    parser.add_argument('--keep_pv',
                        choices=('True', 'False'),
                        help='Whether to keep plausible values',
                        required=False)
    parser.add_argument('--student_id', type=str,
                        help='Name of student id column (common in responses_# and student_exog df)',
                        required=False)
    parser.add_argument('--school_id', type=str,
                        help='name of school id column (in student_exog df)',
                        required=False)
    parser.add_argument('--estimates_in', type=str,
                        help='CSV file with irt estimates for exam 0',
                        required=False)
    parser.add_argument('--responses_in', type=str,
                        help='CSV file with item responses for exam 0',
                        required=False)
    parser.add_argument('--estimates_out', type=str,
                        help='CSV file with irt estimates for exam 1',
                        required=False)
    parser.add_argument('--responses_out', type=str,
                        help='CSV file with item responses for exam 1',
                        required=False)
    parser.add_argument('--student_exog', type=str,
                        help='',
                        required=False)
    parser.add_argument('--fixed_effects', type=str,
                        help='Comma-separated list of names of exog variables (in student_exog file) '
                             'that are to be included in the fixed effects part of latent regression',
                        required=False)
    parser.add_argument('-o','--out_files', type=str,
                        help='Prefix used for CSV files that are saved by pvreg',
                        required=False)

    args = parser.parse_args()

    if args.conf_file:
        with open(args.conf_file) as config_file:
            config = json.load(config_file)

        parser.set_defaults(**config)

    args = parser.parse_args()

    required_list = ['path',
                     'student_id', 'school_id',
                     'estimates_in', 'responses_in',
                     'estimates_out', 'responses_out',
                     'student_exog', 'fixed_effects', 'out_files']
    for action in parser._actions:
        if (action.dest in required_list) and (getattr(args, action.dest) is None):
            action.required = True

    args = parser.parse_args()

    if not args.out_path:
        args.out_path=args.path

    run_pv_reg(args)


if __name__ == "__main__":
    mp.freeze_support() # seems not to change anything - if confirmed delete the line
    main()
