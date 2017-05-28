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
 
