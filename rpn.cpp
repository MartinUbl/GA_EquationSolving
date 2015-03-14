#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "stack.h"
#include "rpn.h"
#include "main.h"

/**
 * Builds element with specified type
 * - there's nothing special in this function, just allocate memory, set type and that's it
 */
rpn_element* rpn_build_element(enum rpn_token_type type)
{
    rpn_element* tmp = (rpn_element*)malloc(sizeof(rpn_element));
    memset(tmp, 0, sizeof(rpn_element));

    tmp->type = type;

    return tmp;
}

/**
 * Applies function (supplied as token identifier) to supplied value and returns it
 */
static double rpn_apply_function(enum supported_functions func, double value)
{
    /* this may fit into some sort of array with function pointers, but I use this
     * because I think, I have better control about what's going on here, and also,
     * as you can see, the argument is not always the same - this may be solved by
     * writing own "intermediate" function which would call the math.h function */

    switch (func)
    {
        case FUNC_ABS:
            return fabs(value);
        case FUNC_EXP:
            return exp(value);
        case FUNC_SIN:
            return sin(value);
        case FUNC_COS:
            return cos(value);
        case FUNC_TAN:
            return tan(value);
        case FUNC_COTAN:
            return 1.0/tan(value);
        case FUNC_ASIN:
            return asin(value);
        case FUNC_ACOS:
            return acos(value);
        case FUNC_ATAN:
            return atan(value);
        case FUNC_ACOTAN:
            return atan(1.0 / value);
        case FUNC_LOG10:
            return log10(value);
        case FUNC_LN:
            return log(value);
        case FUNC_SINH:
            return sinh(value);
        case FUNC_COSH:
            return cosh(value);
        case FUNC_TANH:
            return tanh(value);
        case FUNC_TODEG:
            return value*180.0 / M_PI;
        case FUNC_TORAD:
            return value*M_PI / 180.0;

        default:
        case FUNC_UNSUPPORTED:
            return value;
    }
}

/**
 * Pops two values from stacks and returns them to memory, where supplied pointers point at
 * - also performs cleanup
 */
static void rpn_pop_two_values(c_stack *stck, double *one, double *two)
{
    rpn_element* el;

    el = (rpn_element*)stck_pop(stck);
    *one = el->value.as_double;
    free(el);
    el = (rpn_element*)stck_pop(stck);
    *two = el->value.as_double;
    free(el);
}

/**
 * Helper function for applying operator (supplied as enum type) to N operands on stack
 * the amount of operands is decided in switch body, as well, as their order
 */
static double rpn_apply_operator(enum operator_type op, c_stack *stck)
{
    double one, two;

    switch (op)
    {
        case OP_ADD:
        {
            rpn_pop_two_values(stck, &one, &two);
            return one + two;
        }
        case OP_SUBTRACT:
        {
            rpn_pop_two_values(stck, &one, &two);
            return two - one;
        }
        case OP_MULTIPLY:
        {
            rpn_pop_two_values(stck, &one, &two);
            return one * two;
        }
        case OP_DIVIDE:
        {
            rpn_pop_two_values(stck, &one, &two);
            return two / one;
        }
        case OP_EXP_RAISE:
        {
            rpn_pop_two_values(stck, &one, &two);
            return pow(two, one);
        }
        default:
        {
            return 0.0;
        }
    }
}

/**
 * Helper function for cloning element in stack to be able to destroy the copy
 * instead of original element (we need to reuse the original in more evaluations)
 */
static rpn_element* rpn_clone(rpn_element* source)
{
    rpn_element *tmp = (rpn_element*)malloc(sizeof(rpn_element));
    memcpy(tmp, source, sizeof(rpn_element));
    return tmp;
}

/**
 * Evaluates RPN stack supplued in argument, also considers argument value supplied
 * NOTE: for evaluating using more variables, we would need to assign "IDs" to every
 * variable and supply variable value map with value for every variable. But since we
 * deal with only one variable, passing one variable value is enough
 */
double rpn_evaluate_stack(c_stack* stck, double* variable_values)
{
    int stck_pos;
    double tmp_val;
    rpn_element *rpn_el, *rpn_el_tmp;
    c_stack* rpn_stack;

    /* create new stack used for evaluating */
    rpn_stack = stck_create(stck->curr+1);
    stck_pos = 0;

    /* go through every element in supplied stack */
    while (stck_pos <= stck->curr)
    {
        /* here we access elements in supplied stack as if it's only dynamic array
         * it used to work as stack in preprocessing method (shunting yard algorhitm) */
        rpn_el = (rpn_element*)stck_get(stck, stck_pos);

        /* constants are just pushed to working stack */
        if (rpn_el->type == RPN_TOKEN_CONST)
        {
            stck_push(rpn_stack, rpn_clone(rpn_el));
        }
        /* variables are evaluated using supplied value and pushed to stack */
        else if (rpn_el->type == RPN_TOKEN_VARIABLE)
        {
            rpn_el = rpn_clone(rpn_el);

            rpn_el->type = RPN_TOKEN_CONST;
            rpn_el->value.as_double = variable_values[rpn_el->value.as_variable];

            stck_push(rpn_stack, rpn_el);
        }
        /* processing function will pop value(s) from stack, evaluate them and push result back */
        else if (rpn_el->type == RPN_TOKEN_FUNCTION)
        {
            rpn_el_tmp = (rpn_element*)stck_pop(rpn_stack);
            tmp_val = rpn_apply_function((supported_functions) rpn_el->value.as_function, rpn_el_tmp->value.as_double);
            rpn_el_tmp->value.as_double = tmp_val;

            stck_push(rpn_stack, rpn_el_tmp);
        }
        /* processing operator will pop value(s) from stack, evaluate the expression and push result back */
        else if (rpn_el->type == RPN_TOKEN_OPERATOR)
        {
            tmp_val = rpn_apply_operator((operator_type) rpn_el->value.as_operator, rpn_stack);

            rpn_el_tmp = rpn_build_element(RPN_TOKEN_CONST);
            rpn_el_tmp->value.as_double = tmp_val;

            stck_push(rpn_stack, rpn_el_tmp);
        }

        /* move to next element */
        stck_pos++;
    }

    /* the last element left on stack is our result */
    rpn_el = (rpn_element*)stck_pop(rpn_stack);

    /* take the result and free element */
    tmp_val = (rpn_el != NULL) ? rpn_el->value.as_double : 0.0;
    free(rpn_el);

    /* and clean up the stack we used */
    stck_destroy(rpn_stack);

    return tmp_val;
}
