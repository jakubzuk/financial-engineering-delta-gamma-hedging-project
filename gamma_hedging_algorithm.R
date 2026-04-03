# Funkcja przeprowadzajaca gamma-hedging

# epsilon - kontroluje pasmo (band), transaction_fee - koszt (procent) transakcyjny,
# S0 - poczatkowa wartosc aktywa bazowego, mu - dryf, pricing_sigma - zmiennosc uzyta do wyceny w portfelu replikacyjnym
# market_sigma - zmiennosc wedlug ktorej generujemy trajektorie aktywa bazowego,
# r - stopa procentowa bez ryzyka, t - czas do zapadniecia opcji,
# n - liczba (dat) notowan opcji, m - liczba powtorzen algorytmu
# K - cena wykonania replikowanej opcji, option_premium - cena replikowanej opcji (jesli NA,
# to bierzemy ze wzoru BS), real_prices - notowania aktywa bazowego (jesli NA,
# to generujemy wedlug funkcji gbm_simulation), strike_premium - K*, czyli parametr
# ktory determinuje cene wykonania opcji binarnych ktorymi zabezpieczamy
# beta_upper_limit - maksymalna wartosc bety (ilosci opcji binarnych)

BS_gamma_hedging_algorithm = function(epsilon, transaction_fee, S0, mu, 
                                               pricing_sigma, market_sigma, r, t, 
                                               n, m, K, option_premium = NA, real_prices = NA, 
                                               option_type = 'C', strike_premium = 50, beta_upper_limit = 10e+8) {
  
  gbm_simulation <- function(S0, mu, sigma, t, n, M) {
    dt <- t / n
    
    brownian_increments <- matrix(rnorm(M * n, mean = 0, sd = sqrt(dt)), nrow = M, ncol = n)
    brownian_increments_transformed <- (mu - sigma^2 / 2) * matrix(rep(1:n, M)*dt, nrow = M, byrow = TRUE) + sigma * t(apply(brownian_increments, 1, cumsum))
    brownian_increments_transformed <- cbind(rep(1, M), exp(brownian_increments_transformed))
    
    simulated_price_paths <- S0 * t(brownian_increments_transformed)
    
    return(simulated_price_paths)
  }
  
  BS_price = function(S, K, r, sigma, t_total, t, option_type) {
    d1 = (log(S/K) + (r + sigma^2 / 2)*(t_total - t))/(sigma * (t_total - t)^(1/2))
    d2 = d1 - sigma * (t_total - t)^(1/2)
    
    if (option_type == 'C') {return(S*pnorm(d1) - K*exp(-r*(t_total - t))*pnorm(d2))}
    else {return(-S*pnorm(-d1) + K*exp(-r*(t_total - t))*pnorm(-d2))}
  }
  
  BS_delta = function(S, K, r, sigma, t_total, t, option_type) {
    if (option_type == 'C') {return(pnorm((log(S/K) + (r + sigma^2 / 2)*(t_total - t))/(sigma * (t_total - t)^(1/2))))}
    else {return(pnorm((log(S/K) + (r + sigma^2 / 2)*(t_total - t))/(sigma * (t_total - t)^(1/2))) - 1)}
  }
  
  BS_gamma = function(S, K, r, sigma, t_total, t) {
    d1 = (log(S/K) + (r + sigma^2 / 2)*(t_total - t))/(sigma * (t_total - t)^(1/2))
    return(dnorm(d1)/(sigma * S * sqrt(t_total - t)))
  }
  
  BS_binary_price = function(S, K, r, sigma, t_total, t, option_type) {
    d1 = (log(S/K) + (r + sigma^2 / 2)*(t_total - t))/(sigma * (t_total - t)^(1/2))
    d2 = d1 - sigma * (t_total - t)^(1/2)
    
    if (option_type == 'C') {return(exp(-r * (t_total - t))*pnorm(d2))}
    else {return(exp(-r * (t_total - t))*(1 - pnorm(d2)))}
  }
  
  BS_binary_delta = function(S, K, r, sigma, t_total, t, option_type) {
    d1 = (log(S/K) + (r + sigma^2 / 2)*(t_total - t))/(sigma * (t_total - t)^(1/2))
    d2 = d1 - sigma * (t_total - t)^(1/2)
    
    if (option_type == 'C') {return(exp(-r*(t_total - t)) * dnorm(d2) / (sigma * S * (t_total - t)^(1/2)))}
    else {return(-exp(-r * (t_total - t)) * dnorm(d2) / (sigma * S * (t_total - t)^(1/2)))}
  }
  
  BS_binary_gamma = function(S, K, r, sigma, t_total, t, option_type) {
    d1 = (log(S/K) + (r + sigma^2 / 2)*(t_total - t))/(sigma * (t_total - t)^(1/2))
    d2 = d1 - sigma * (t_total - t)^(1/2)
    
    if (option_type == 'C') {return( - exp(-r*(t_total - t)) * d1 * dnorm(d2) / (sigma^2 * S^2 * (t_total - t)))}
    else {return(exp(-r*(t_total - t)) * d1 * dnorm(d2) / (sigma^2 * S^2 * (t_total - t)))}
  }
  
  if ((n < 2) | (prod(is.na(real_prices)) & (n > t*250)) | ((!prod(is.na(real_prices))) & (n > t*length(real_prices)))) {
    message('n (częstość aktualizowania portfela) za niska lub za wysoka')
    stop()
  }
  
  alphas_matrix = matrix(nrow = n, ncol = m)
  betas_matrix = matrix(nrow = n, ncol = m)
  cashs_matrix = matrix(nrow = n, ncol = m)
  option_values_matrix = matrix(nrow = n, ncol = m)
  final_price_diffs = numeric(m)
  transactions = numeric(m)
  transaction_costs = numeric(m)
  
  for (j in 1:m) {
    
    if (!!prod(is.na(real_prices))) {price_path = c(gbm_simulation(S0, mu, market_sigma, t, 250, 1))}
    else {price_path = real_prices}

    tstep = 0
    dt = t/n
    
    if (is.na(option_premium)) {
      cashs_matrix[1, j] = BS_price(price_path[1], K, r, pricing_sigma, t, tstep, option_type)
      option_values_matrix[1, j] = BS_price(price_path[1], K, r, pricing_sigma, t, tstep, option_type)
    }
    else {
      cashs_matrix[1, j] = option_premium
      option_values_matrix[1, j] = option_premium
    }
    
    delta_B = BS_binary_delta(price_path[1], price_path[1] + strike_premium, r, pricing_sigma, t, tstep, 'C') + BS_binary_delta(price_path[1], price_path[1] - strike_premium, r, pricing_sigma, t, tstep, 'P')
    beta = min(beta_upper_limit, BS_gamma(price_path[1], K, r, pricing_sigma, t, tstep) / (BS_binary_gamma(price_path[1], price_path[1] + strike_premium, r, pricing_sigma, t, tstep, 'C') + BS_binary_gamma(price_path[1], price_path[1] - strike_premium, r, pricing_sigma, t, tstep, 'P')))
    
    betas_matrix[1, j] = beta
    alphas_matrix[1, j] = BS_delta(price_path[1], K, r, pricing_sigma, t, tstep, option_type) - delta_B * beta
    cashs_matrix[1, j] = cashs_matrix[1, j] - alphas_matrix[1, j] * price_path[1] - beta * (BS_binary_price(price_path[1], price_path[1] + strike_premium, r, pricing_sigma, t, tstep, 'C') + BS_binary_price(price_path[1], price_path[1] - strike_premium, r, pricing_sigma, t, tstep, 'P'))
    
    latest_binary_call_strike = price_path[1] + strike_premium
    latest_binary_put_strike = price_path[1] - strike_premium
    
    transactions[j] = transactions[j] + 1
    transaction_costs[j] = transaction_costs[j] + transaction_fee * abs(alphas_matrix[1, j]) * price_path[1] + transaction_fee * abs(betas_matrix[1, j]) * (BS_binary_price(price_path[1], price_path[1] + strike_premium, r, pricing_sigma, t, tstep, 'C') + BS_binary_price(price_path[1], price_path[1] - strike_premium, r, pricing_sigma, t, tstep, 'P'))
    
    for (i in 1:(n - 1)) {
      
      tstep = tstep + dt
      
      option_values_matrix[i + 1, j] = BS_price(price_path[floor(length(price_path)*tstep) + 1], K, r, pricing_sigma, t, tstep, option_type)
      
      new_delta_B = BS_binary_delta(price_path[floor(length(price_path)*tstep) + 1], price_path[floor(length(price_path)*tstep) + 1] + strike_premium, r, pricing_sigma, t, tstep, 'C') + BS_binary_delta(price_path[floor(length(price_path)*tstep) + 1], price_path[floor(length(price_path)*tstep) + 1] - strike_premium, r, pricing_sigma, t, tstep, 'P')
      new_beta = min(beta_upper_limit, BS_gamma(price_path[floor(length(price_path)*tstep) + 1], K, r, pricing_sigma, t, tstep) / (BS_binary_gamma(price_path[floor(length(price_path)*tstep) + 1], price_path[floor(length(price_path)*tstep) + 1] + strike_premium, r, pricing_sigma, t, tstep, 'C') + BS_binary_gamma(price_path[floor(length(price_path)*tstep) + 1], price_path[floor(length(price_path)*tstep) + 1] - strike_premium, r, pricing_sigma, t, tstep, 'P')))
      new_alpha = BS_delta(price_path[floor(length(price_path)*tstep) + 1], K, r, pricing_sigma, t, tstep, option_type) - new_delta_B * new_beta
      
      if (abs(new_alpha - alphas_matrix[i, j]) < epsilon) {
        alphas_matrix[i + 1, j] = alphas_matrix[i, j]
        betas_matrix[i + 1, j] = betas_matrix[i, j]
        cashs_matrix[i + 1, j] = cashs_matrix[i, j] * exp(r*dt)
        next
      }
      
      alphas_matrix[i + 1, j] = new_alpha
      betas_matrix[i + 1, j] = new_beta
      
      new_binaries_price = BS_binary_price(price_path[floor(length(price_path)*tstep) + 1], price_path[floor(length(price_path)*tstep) + 1] + strike_premium, r, pricing_sigma, t, tstep, 'C') + BS_binary_price(price_path[floor(length(price_path)*tstep) + 1], price_path[floor(length(price_path)*tstep) + 1] - strike_premium, r, pricing_sigma, t, tstep, 'P')
      old_binaries_price = BS_binary_price(price_path[floor(length(price_path)*tstep) + 1], latest_binary_call_strike, r, pricing_sigma, t, tstep, 'C') + BS_binary_price(price_path[floor(length(price_path)*tstep) + 1], latest_binary_put_strike, r, pricing_sigma, t, tstep, 'P')
      cashs_matrix[i + 1, j] = cashs_matrix[i, j] * exp(r*dt) - (alphas_matrix[i + 1, j] - alphas_matrix[i, j])*price_path[floor(length(price_path)*tstep) + 1] + betas_matrix[i, j] * old_binaries_price - betas_matrix[i + 1, j] * new_binaries_price
      
      latest_binary_call_strike = price_path[floor(length(price_path)*tstep) + 1] + strike_premium
      latest_binary_put_strike = price_path[floor(length(price_path)*tstep) + 1] - strike_premium
      
      transactions[j] = transactions[j] + 1
      transaction_costs[j] = transaction_costs[j] + transaction_fee * abs(alphas_matrix[i, j] - alphas_matrix[i + 1, j]) * price_path[floor(length(price_path)*tstep) + 1] + transaction_fee * abs(betas_matrix[i, j]) * old_binaries_price + transaction_fee * abs(betas_matrix[i + 1, j])* new_binaries_price
    }
    

    if (option_type == 'C') {payoff = max(price_path[length(price_path)] - K, 0)}
    else {payoff = max(K - price_path[length(price_path)], 0)}
    
    binary_call_payoff = ifelse(price_path[length(price_path)] > latest_binary_call_strike, 1, 0)
    binary_put_payoff = ifelse(price_path[length(price_path)] < latest_binary_put_strike, 1, 0)
    
    portfolio = cashs_matrix[n, j] * exp(r*(t - tstep)) + alphas_matrix[n, j] * price_path[length(price_path)] + betas_matrix[n, j] * (binary_call_payoff + binary_put_payoff)
    

    final_price_diffs[j] = portfolio - payoff
  }
  
  return(list(final_price_diffs, alphas_matrix, cashs_matrix, transactions, transaction_costs, option_values_matrix, betas_matrix)) 
}


# Epsilony jakie rozwazamy (dla innych, tj. wiekszych od 0.11 slabe rezultaty)
epsilon_range = c(0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11)
# K* jakie rozwazamy (dla wiekszych problem, bo zdaza sie ze cena spada ponizej 700, i mamy negatywny strike dla binarnych putow, poza tym
# wystarczajaco dobre wyniki zachodza juz dla nizszych)
strike_premium_range = c(50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700)

# Parametry
transaction_fee = 0.004
S0 = 2270
mu = 0.1
pricing_sigma = 0.2
market_sigma = 0.2
r = 0.05
t = 1
n = 250
M = 100
K = 1800
beta_limit = 10^4

# Lista z rezultatami gamma hedgingu dla kombinacji epsilonow i K*
gamma_hedging_data_list = list()
# Macierze ze statystykami opisujacymi wyniki dla kombinacji epsilonow i K*
gamma_hedging_means_matrix = matrix(nrow = length(epsilon_range), ncol = length(strike_premium_range))
gamma_hedging_sds_matrix = matrix(nrow = length(epsilon_range), ncol = length(strike_premium_range))
gamma_hedging_mean_betas_matrix = matrix(nrow = length(epsilon_range), ncol = length(strike_premium_range))

rownames(gamma_hedging_means_matrix) = epsilon_range
colnames(gamma_hedging_means_matrix) = strike_premium_range
rownames(gamma_hedging_sds_matrix) = epsilon_range
colnames(gamma_hedging_sds_matrix) = strike_premium_range
rownames(gamma_hedging_mean_betas_matrix) = epsilon_range
colnames(gamma_hedging_mean_betas_matrix) = strike_premium_range

for (i in 1:length(epsilon_range)) {
  gamma_hedging_data_list[[i]] = list()
}

# Ta zmienna sluzy do tego, ze czasami (przy bardzo rzadkich trajektoriach) w gamma-hedgingu, przez zabezpieczanie opcjami binarnymi
# mozemy uzyskac duzy profit, ktory zaburza statystyki (sztucznie podnosi srednia roznice miedzy portfelem a payoffem opcji), wiec
# przyjmujac np. positive_limit = 100 mozemy wyeliminowac takie obserwacje (dla positive_limit = Inf rozwazamy wszystkie obserwacje)
positive_limit = Inf

# Wyznaczanie danych (to sie chwilke kompiluje, u mnie dla M = 100 tak z 5 minut)
for (i in 1:length(epsilon_range)) {
  for (j in 1:length(strike_premium_range)) {
    gamma_hedging_data_list[[i]][[j]] = BS_gamma_hedging_algorithm(epsilon_range[i], transaction_fee, S0,
                                                                   mu, pricing_sigma, market_sigma, r, t,
                                                                   n, M, K, strike_premium = strike_premium_range[j],
                                                                   beta_upper_limit = beta_limit)
    
    gamma_hedging_means_matrix[i, j] = mean((gamma_hedging_data_list[[i]][[j]][[1]] - gamma_hedging_data_list[[i]][[j]][[5]])[gamma_hedging_data_list[[i]][[j]][[1]] - gamma_hedging_data_list[[i]][[j]][[5]] < positive_limit])
    gamma_hedging_sds_matrix[i, j] = sd((gamma_hedging_data_list[[i]][[j]][[1]] - gamma_hedging_data_list[[i]][[j]][[5]])[gamma_hedging_data_list[[i]][[j]][[1]] - gamma_hedging_data_list[[i]][[j]][[5]] < positive_limit])
    gamma_hedging_mean_betas_matrix[i, j] = mean(gamma_hedging_data_list[[i]][[j]][[7]][, gamma_hedging_data_list[[i]][[j]][[1]] - gamma_hedging_data_list[[i]][[j]][[5]] < positive_limit])
  }
}


library(ggplot2)
library(reshape2)

df_sds = melt(t(gamma_hedging_sds_matrix))
df_sds$Var1 = factor(df_sds$Var1)
df_sds$Var2 = factor(df_sds$Var2)

# Mapa ciepla odchylen standardowych roznic miedzy wartoscia portfela a payoffem opcji
sds_plot = ggplot(df_sds, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 1)), size = 3) +  # <-- TU
  coord_fixed() + ylab(bquote(epsilon)) + xlab('K*') + 
  labs(title = "Odchylenie standardowe roznicy miedzy portfelem zabezpieczajacym \na payoffem opcji zabezpieczanej (call, K = 1800), \nGamma-Hedging z danym epsilonem (kontrolujacym band) \ni parametrem K* modyfikujacym strike'i opcji binarnych",
       subtitle = bquote("Parametry rynku: " ~ sigma ~ '= 0.2' ~ mu ~ '= 0.1' ~ 'r = 0.05')) +
  scale_fill_gradientn(
    colors = c("green", "yellow", "orange", "red"),
    values = scales::rescale(c(0, 10, 50, 300))
  ) + scale_y_discrete(limits = rev(levels(df$Var2))) + theme_minimal() + coord_fixed(ratio = 0.75) +
  theme(legend.position = 'none', plot.title = element_text(size = 12), plot.subtitle = element_text(size = 10))
  

df_means = melt(t(gamma_hedging_means_matrix))
df_means$Var1 = factor(df_means$Var1)
df_means$Var2 = factor(df_means$Var2)

# Mapa ciepla srednich roznic miedzy wartoscia portfela a payoffem opcji
means_plot = ggplot(df_means, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 1)), size = 3) +  # <-- TU
  coord_fixed() + ylab(bquote(epsilon)) + xlab('K*') + 
  labs(title = "Srednia roznica miedzy portfelem zabezpieczajacym \na payoffem opcji zabezpieczanej (call, K = 1800), \nGamma-Hedging z danym epsilonem (kontrolujacym band) \ni parametrem K* modyfikujacym strike'i opcji binarnych",
       subtitle = bquote("Parametry rynku: " ~ sigma ~ '= 0.2' ~ mu ~ '= 0.1' ~ 'r = 0.05')) +
  scale_fill_gradientn(
    colors = c("red", "orange", "yellow", "green"),
    values = scales::rescale(c(-500, -75, -30, -10))
  ) + scale_y_discrete(limits = rev(levels(df$Var2))) + theme_minimal() + coord_fixed(ratio = 0.75) +
  theme(legend.position = 'none', plot.title = element_text(size = 12), plot.subtitle = element_text(size = 10))

# Wyswietlenie wykresow
means_plot | sds_plot




# Jesli przyjmiemy epsilon = 0 i beta_upper_limit = 0, to funkcja przeprowadza
# standardowy delta-hedging z codziennym rehedgingiem
delta_hedging_case = BS_gamma_hedging_algorithm(0, transaction_fee, S0, mu, pricing_sigma, market_sigma,
                                                r, t, n, 1000, K, strike_premium = 10, beta_upper_limit = 0)
delta_hedging_cost_included_mean = mean(delta_hedging_case[[1]] - delta_hedging_case[[5]])
delta_hedging_cost_included_sd = sd(delta_hedging_case[[1]] - delta_hedging_case[[5]])
# Srednia wynosi okolo -21.80, a odchylenie standardowe 13.70 (oba z uwzglednionymi kosztami)


