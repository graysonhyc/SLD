@mfunction("_lambda, rs, ybar")
def training(freq=None, amount=None):

    #learning rate
    eta1 = 0.06
    eta2 = 0.001
    eta3 = 1.5

    _lambda = freq(1)# estim of lambda

    month = (mslice[0:58])

    for i in mslice[2:30]:
        e1 = freq(i) - _lambda
        _lambda = _lambda + eta1 * e1
        error_l(i - 1).lvalue = e1 / freq(i)
        end

        for i in mslice[1:60]:
            n = 0
            amount_A = 0
            for j in mslice[1:27]:
                if amount(j, i) != 0:
                    amount_A(n + 1).lvalue = amount(j, i)
                    n = n + 1
                    end
                    end

                    if i == 1:
                        ybar = mean(amount_A)                    # sample mean
                        rs = var(amount_A)                    # sample variance = estim of var parameter
                        #xi = exp(ybar+rs^2/2) ;% mean parameter of log normal dist
                        xi = exp(ybar + rs ** 2 / 2)
                    else:
                        rsnew = var(amount_A)
                        ybarnew = mean(amount_A)
                        xinew = exp(ybarnew + rsnew ** 2 / 2)
                        e2 = rsnew - rs
                        rs = rs + eta2 * e2
                        e3 = xinew - xi
                        xi = xinew + eta3 * e3
                        error_x(i - 1).lvalue = abs((xinew - xi) / xi)
                        end
                        end


                        #figure
                        #1plot(month,abs(error_l))
                        #figure
                        #plot(month,y)
                        #figure
                        plot(month, error_x)
                        end
 
hold(mstring('on'))
t = datetime(2017, 3, 28) + caldays(mslice[1:120])

for k in mslice[1:10]:

    n = 20# number of samples
    lam = freq# Poisson lambda
    mu = xi # lognormal mu parameter
    sig = 0# lognormal sigma parameter
    fn = lambda m: poissrnd(lam, m, 1)# fn(m) returns a Poisson sample of size m
    fx = lambda m: exp(mu + sig * randn(m, 1))# fx(m) returns a lognormal sample of size m
    y1 = zeros(n, 1)
    y3 = zeros(n, 1)
    y5 = zeros(n, 1)
    y7 = zeros(n, 1)

    for i in mslice[1:numel(y1)]:
        nr = fn(1)
        for j in mslice[1:nr]:
            y1(i).lvalue = y1(i) + fx(1)
            y3(i).lvalue = y3(i) + fx(1)
            y5(i).lvalue = y5(i) + fx(1)
            y7(i).lvalue = y7(i) + fx(1)
            end
            end
            n = 10
            lam = freq
            mu = xi       # lognormal mu parameter
            sig = 0.5
            fn = lambda m: poissrnd(lam, m, 1)        # fn(m) returns a Poisson sample of size m
            fx = lambda m: exp(mu + sig * randn(m, 1))        # fx(m) returns a lognormal sample of size m
            y2 = zeros(n, 1)
            y4 = zeros(n, 1)
            y6 = zeros(n, 1)
            y8 = zeros(n, 1)

            for i in mslice[1:numel(y2)]:
                nr = fn(1)
                for j in mslice[1:nr]:
                    y2(i).lvalue = y2(i) + fx(1)
                    y4(i).lvalue = y4(i) + fx(1)
                    y6(i).lvalue = y6(i) + fx(1)
                    y8(i).lvalue = y8(i) + fx(1)

                    end
                    end

                    y(mslice[1:20]).lvalue = y1
                    y(mslice[21:30]).lvalue = y2
                    y(mslice[31:50]).lvalue = y3
                    y(mslice[51:60]).lvalue = y4
                    y(mslice[61:80]).lvalue = y5
                    y(mslice[81:90]).lvalue = y6
                    y(mslice[91:110]).lvalue = y7
                    y(mslice[111:120]).lvalue = y8

                    acct(mslice[1:30]).lvalue = 700
                    acct(mslice[31:60]).lvalue = 1400
                    acct(mslice[61:90]).lvalue = 2100
                    acct(mslice[91:120]).lvalue = 2800

                    for i in mslice[1:120]:
                        acct(i).lvalue = acct(i) - sum(y(mslice[1:i]))
                        end
                        acct = 20 * acct
                        plot(t, acct)

                        end
 

# function [P, Nl] = mlmc(Lmin,Lmax,N0,eps,mlmc_l, alpha,beta,gamma)
#
# multi-level Monte Carlo estimation
#
# P     = value
# Nl    = number of samples at each level
#
# Lmin  = minimum level of refinement       >= 2
# Lmax  = maximum level of refinement       >= Lmin
# N0    = initial number of samples         > 0
# eps   = desired accuracy (rms error)      > 0 
#
# alpha -> weak error is  O(2^{-alpha*l})
# beta  -> variance is    O(2^{-beta*l})
# gamma -> sample cost is O(2^{gamma*l})    > 0
#
# if alpha, beta are not positive then they will be estimated
#
# mlmc_l = function for level l estimator 
#
# sums = mlmc_fn(l,N)     low-level routine
#
# inputs:  l = level
#          N = number of paths
#
# output: sums(1) = sum(Y)
#         sums(2) = sum(Y.^2)
#         where Y are iid samples with expected value:
#         E[P_0]           on level 0
#         E[P_l - P_{l-1}] on level l>0

@mfunction("P, Nl")
def mlmc(Lmin=None, Lmax=None, N0=None, eps=None, mlmc_l=None, alpha_0=None, beta_0=None, gamma=None):

    #
    # check input parameters
    #
    if (Lmin < 2):
        error(mstring('error: needs Lmin >= 2'))
        end

        if (Lmax < Lmin):
            error(mstring('error: needs Lmax >= Lmin'))
            end

            if (N0 <= 0 or eps <= 0 or gamma <= 0):
                error(mstring('error: needs N>0, eps>0, gamma>0 \\n'))
                end

                #
                # initialisation
                #
                alpha = max(0, alpha_0)
                beta = max(0, beta_0)

                theta = 0.25

                L = Lmin

                Nl(mslice[1:L + 1]).lvalue = 0
                suml(mslice[1:2], mslice[1:L + 1]).lvalue = 0
                dNl(mslice[1:L + 1]).lvalue = N0

                while sum(dNl) > 0:

                    #
                    # update sample sums
                    #
                    for l in mslice[0:L]:
                        if dNl(l + 1) > 0:
                            sums = feval(mlmc_l, l, dNl(l + 1))
                            Nl(l + 1).lvalue = Nl(l + 1) + dNl(l + 1)
                            suml(1, l + 1).lvalue = suml(1, l + 1) + sums(1)
                            suml(2, l + 1).lvalue = suml(2, l + 1) + sums(2)
                            end
                            end

                            #
                            # compute absolute average and variance
                            #
                            ml = abs(suml(1, mslice[:]) /eldiv/ Nl)
                            Vl = max(0, suml(2, mslice[:]) /eldiv/ Nl - ml **elpow** 2)

                            #
                            # fix to cope with possible zero values for ml and Vl
                            # (can happen in some applications when there are few samples)
                            #
                            for l in mslice[3:L + 1]:
                                ml(l).lvalue = max(ml(l), 0.5 * ml(l - 1) / 2 ** alpha)
                                Vl(l).lvalue = max(Vl(l), 0.5 * Vl(l - 1) / 2 ** beta)
                                end

                                #
                                # use linear regression to estimate alpha, beta if not given
                                #
                                if alpha_0 <= 0:
                                    A = repmat((mslice[1:L]).cT, 1, 2) **elpow** repmat(mslice[1:-1:0], L, 1)
                                    x = A; print x
                                    log2(ml(mslice[2:end])).cT

                                    alpha = max(0.5, -x(1))
                                    end

                                    if beta_0 <= 0:
                                        A = repmat((mslice[1:L]).cT, 1, 2) **elpow** repmat(mslice[1:-1:0], L, 1)
                                        x = A; print x
                                        log2(Vl(mslice[2:end])).cT

                                        beta = max(0.5, -x(1))
                                        end
                                        #
                                        # set optimal number of additional samples
                                        #
                                        Cl = 2. ** (gamma * (mslice[0:L]))
                                        Ns = ceil(sqrt(Vl /eldiv/ Cl) * sum(sqrt(Vl *elmul* Cl)) / ((1 - theta) * eps ** 2))
                                        dNl = max(0, Ns - Nl)
                                        #
                                        # if (almost) converged, estimate remaining error and decide 
                                        # whether a new level is required
                                        #
                                        if sum(dNl > 0.01 * Nl) == 0:
                                            rem = ml(L + 1) / (2 ** alpha - 1)

                                            if rem > sqrt(theta) * eps:
                                                if (L == Lmax):
                                                    fprintf(1, mstring('*** failed to achieve weak convergence *** \\n'))
                                                else:
                                                    L = L + 1
                                                    Vl(L + 1).lvalue = Vl(L) / 2 ** beta
                                                    Nl(L + 1).lvalue = 0
                                                    suml(mslice[1:4], L + 1).lvalue = 0

                                                    Cl = 2. ** (gamma * (mslice[0:L]))
                                                    Ns = ceil(sqrt(Vl /eldiv/ Cl) * sum(sqrt(Vl *elmul* Cl)) / ((1 - theta) * eps ** 2))
                                                    dNl = max(0, Ns - Nl)
                                                    end
                                                    end
                                                    end
                                                    end

                                                    #
                                                    # finally, evaluate multilevel estimator
                                                    #
                                                    P = sum(suml(1, mslice[:]) /eldiv/ Nl)
                                                    end
 
