import numpy as np
import sympy as sp
import inspect
import types



# Solution Class
class NM_Solution:
    def __init__(self , method , root , max_iteration):
        self.method = method
        self.root = root
        self.max_iteration = max_iteration


# Iteration Limit        
iteration_exit = 100

def error(new , old):
    return np.abs( (new - old) / new ) * 100

# Check for lambda or sympy functions
def check_func(f , sym = False , lamb = False):

    if sym:
        if isinstance(f, sp.Basic):
            return True
        else:
            return False
    if lamb:
        if isinstance(f, types.LambdaType):
            if f.__name__ == "<lambda>":
                return True
            else:
                return False
        else:
            return False

    if isinstance(f, types.LambdaType) or isinstance(f, sp.Basic):
        return True
    else:
        return False

# Lambda Function to Symbolic Function (sympy)
def lambda2sym(lambda_function , symble = None):
    if not check_func(lambda_function , lamb = True):
        raise ValueError("Not a lambda function.")
    
    if not symble:
        symble = 'x'

    if check_func(symble , sym=True):
        pass
    else:
        # Define the Symble for Sympy function
        exec(f"{symble}" + " = sp.symbols('" + f"{symble}" + "')")
    
    # Get the source code of the lambda function as a string
    lambda_function_string = inspect.getsource(lambda_function)

    # Split the string using ':' as the delimiter and get the second part
    expression_part = lambda_function_string.split(':')[1].strip()

    # Return as Sympy function
    # return eval(expression_part)
    return sp.sympify(expression_part)


#########################################################
#                                                       #
#   Methods                                             #
#                                                       #
#########################################################

# Bisection Method

def bisection_method( func , x_lower , x_upper , true_value = None , terminal_error = 5 , iteration_limit=None , significant_digit = None , analysis = False):
    # iteration_exit = 100

    if check_func(func , sym=True):
        raise ValueError("Symbolic Function is not supported yet in <" + __package__ + "> .")

    # Creating terminal_error for significant_digit
    if significant_digit:
        terminal_error = 0.5 * ( 10**(2 - significant_digit) )

    if func(x_lower)*func(x_upper) < 0:
        relative_error = 100
        x_r_old = 0
        iteration = 1
        while relative_error > terminal_error:
            x_r = (x_lower + x_upper) / 2
            if func(x_lower)*func(x_r) < 0:
                x_upper = x_r
            elif func(x_r)*func(x_upper) < 0:
                x_lower = x_r
            elif func(x_r) == 0:
                break
            
            relative_error = error( x_r , x_r_old )
            x_r_old = x_r
            iteration += 1

            if not iteration_limit:
                if iteration > iteration_exit:
                    print("Exceeding <iteration_limit>. Increase the <iteration_limit> to do longer itertation ")
                    break 
            elif iteration > iteration_limit:
                break
                      
        
    else:
        raise ValueError("Guess is not current.")
    
    if analysis:
        return NM_Solution( "Bisection Method" , x_r , iteration )
    else:
        return x_r

# False-Position Method

def false_position_method( func , x_lower , x_upper , true_value = None , terminal_error = 5 , iteration_limit=None , significant_digit = None , analysis = False):
    # iteration_exit = 100

    if check_func(func , sym=True):
        raise ValueError("Symbolic Function is not supported yet in <" + __package__ + "> .")

    # Creating terminal_error for significant_digit
    if significant_digit:
        terminal_error = 0.5 * ( 10**(2 - significant_digit) )

    if func(x_lower)*func(x_upper) < 0:
        relative_error = 100
        x_r_old = 0
        iteration = 1
        while relative_error > terminal_error:

            x_r = x_upper - ( func(x_upper) * ( x_upper - x_lower ) / ( func(x_upper) - func(x_lower) ) ) 
            
            if func(x_lower)*func(x_r) < 0:
                x_upper = x_r
            elif func(x_r)*func(x_upper) < 0:
                x_lower = x_r
            elif func(x_r) == 0:
                break
            
            relative_error = error(x_r , x_r_old)
            x_r_old = x_r
            iteration += 1

            if not iteration_limit:
                if iteration > iteration_exit:
                    print("Exceeding <iteration_limit>. Increase the <iteration_limit> to do longer itertation ")
                    break 
            elif iteration > iteration_limit:
                break
                      
        
    else:
        raise ValueError("Guess is not current.")
    
    if analysis:
        return NM_Solution( "False Position Method" , x_r , iteration )
    else:
        return x_r



# Fixed-Point Iteration

def fixed_point_iteration( func , variable , x_guess , gofx_format = True , true_value = None , terminal_error = 5 , iteration_limit=None , significant_digit = None , analysis = False ):
    
    if not check_func( variable , sym=True ):
        variable = sp.symbols(variable)

    if not check_func( func , sym=True ):
        func = lambda2sym( func , variable )

    if not gofx_format:
        func = func + variable

    # Creating terminal_error for significant_digit
    if significant_digit:
        terminal_error = 0.5 * ( 10**(2 - significant_digit) )

    func_dif = sp.diff(func)

    relative_error = 100
    x_old = x_guess
    iteration = 1
    while relative_error > terminal_error:
        
        if func_dif.subs(variable , x_old) < 1:
            x_new = func.subs(variable , x_old)
            relative_error = error(x_new , x_old)

            x_old = x_new
            iteration += 1

            if not iteration_limit:
                if iteration > iteration_exit:
                    print("Exceeding <iteration_limit>. Increase the <iteration_limit> to do longer itertation ")
                    break 
            elif iteration > iteration_limit:
                break
                      
        
        else:
            print("<Error:> Function will not converge.")
            break
        
    if analysis:
        return NM_Solution( "Fixed Point Iteration" , x_old , iteration )
    else:
        return x_new

    
        pass


# Newton-Raphson