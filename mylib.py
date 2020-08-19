#comment added to explore pull in git lesson

def is_prime(num):
    if num == 2:
        return True
    if num < 2 or num % 2 == 0:
        return False
    for divisor in range(3, num//2+1, 2):
        if num % divisor == 0:
            return False
    return True
####################################
def find_divisors(n):
    lst = []
    for num in range(1, n + 1):
        if n % num == 0:
            lst.append(num)
    return lst

##################################
def delete_listdupes(lst):   #This is the often referred to simple set trick to dedup a list,
    new_lst = list(set(lst))  #string, other simple collections
    return new_lst

x = [1, 3, 1, 1, 2, 2, 1, 3, 3, 1, 3, 1, 3, 3, 3, 1, 3, 3, 3, 1]

###################################
def delete_dups(d):  #for deduping a dictionary
    dcopy = d.copy()
    result = {}
        
    for k, v in dcopy.items():
        if v not in result.values():
            result[k] = v
    return result

##################################
def mean(lst):
    return sum(lst) / len(lst)

################################
def median(lst):
    lst_sorted = sorted(lst)
    
    if len(lst) % 2 == 1:
        median_idx = int(len(lst) / 2)
        return lst_sorted[median_idx]
    else:
        higher_mid = lst_sorted[int(len(lst)/2)]
        lower_mid = lst_sorted[int(len(lst)/2)-1]

        return (higher_mid + lower_mid) / 2


######################################
def mode(lst):
    most_occurring = lst[0]

    for item in lst[1:]:
        if lst.count(item) > lst.count(most_occurring):
            most_occurring = item

    return most_occurring   
#####################################
def five_summary(lst):
    # Sort list, find median
    sorted_lst = sorted(lst)
    med = median(sorted_lst)

    # If odd length, find index of median, partition data
    if len(lst) % 2 != 0:
        med_idx = (len(lst) / 2 + .5) - 1
        low_subset = sorted_lst[:int(med_idx)]
        high_subset = sorted_lst[int(med_idx + 1):]

    # If even length, find index of median, partition data
    else:
        idx1 = int(len(lst) / 2 - 1)
        idx2 = int(idx1 + 1)

        low_subset = sorted_lst[0:(idx1 + 1)]
        high_subset = sorted_lst[idx2:]

    # Define Q1, Q3
    q1 = median(low_subset)
    q3 = median(high_subset)

    # Return five number summary in a tuple
    return min(sorted_lst), q1, med, q3, max(sorted_lst)
######################################################
def variance(lst):  #population
    """Returns a float which represents the variance of a given list, this
    function does not implement Bessel's correction.

    Parameters
    ----------
    lst : list
        A list which represents a dataset that is comprised entirely of numerical values

    Returns
    -------
    variance : float
        A float which represents the variance of the dataset
    """
    the_mean = mean(lst)
    v_sum = 0

    # For each value in the list, find the squared difference
    for value in lst:
        v_sum += (value - the_mean)**2

    # Multiply by 1 / n
    return v_sum / len(lst)

 ############################################
def variance_sample(lst):  #sample only, reduce to n-1 for divisor
    the_mean = mean(lst)
    v_sum = 0

    # For each value in the list, find the squared difference
    for value in lst:
        v_sum += (value - the_mean)**2

    # Multiply by 1 / n-1
    return v_sum / (len(lst)-1)
##############################################################
def factorial(n):
    prod = 1
    for num in range(1, n+1):
        prod *= num
    return prod
###############################################################
def permutations(n, r):
    return factorial(n) / factorial(n-r)

######################################
def combinations(n, k): 
    return int(factorial(n) / (factorial(n-k) * factorial(k)))
print(combinations(15,2))
######################################
def binomial_pmf(n, k, p=0.5):
    return combinations(n, k) * (p**k) * (1-p)**(n-k)

######################################
def binomial_cdf(n, k_high, p=0.5):
    cumulative = 0

    for k in range(k_high+1):
        cumulative += binomial_pmf(n, k, p)
    
    return cumulative

####################################
import math
# Must import math library for this to work
def poisson_dist(occurrence, rate):
    pd = math.e**(-rate) * rate**occurrence / math.factorial(occurrence)
    
    return pd

######################################
