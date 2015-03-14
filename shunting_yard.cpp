#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "stack.h"
#include "rpn.h"
#include "main.h"
#include "shunting_yard.h"

/* helpful macro for routine used when hit some syntax error */
#define ERROR_ROUTINE(a,b) *error=a; *error_ptr=b; stck_destroy(op_stack); stck_destroy(rpn_stack);

static func_match_template func_match[] = {
        { FUNC_ABS, "abs" },
        { FUNC_EXP, "exp" },

        { FUNC_SIN, "sin" },
        { FUNC_COS, "cos" },
        { FUNC_TAN, "tan" },
        { FUNC_COTAN, "cotan" },

        { FUNC_ASIN, "asin" },
        { FUNC_ACOS, "acos" },
        { FUNC_ATAN, "atan" },
        { FUNC_ACOTAN, "acotan" },

        { FUNC_LOG10, "log" },
        { FUNC_LN, "ln" },

        { FUNC_SINH, "sinh" },
        { FUNC_COSH, "cosh" },
        { FUNC_TANH, "tanh" },

        { FUNC_TODEG, "todeg" },
        { FUNC_TORAD, "torad" }
};

/**
 * Helper function to retrieve numeric representation of character digit
 * if not a number, return -1
 */
static int sy_get_number(char chr)
{
    if (chr >= '0' && chr <= '9')
        return chr - '0';

    return -1;
}

/**
 * Helper function to retrieve operator token identifier
 * if not an operator, return -1 (OP_NONE constant)
 */
static int sy_get_operator(char chr)
{
    switch (chr)
    {
        case '+':
            return OP_ADD;
        case '-':
            return OP_SUBTRACT;
        case '*':
            return OP_MULTIPLY;
        case '/':
            return OP_DIVIDE;
        case '^':
            return OP_EXP_RAISE;
    }

    return OP_NONE;
}

/**
 * Helper function to retrieve parenthesis token identifier
 * if not a parenthesis, return -1 (PARENTHESIS_NONE constant)
 */
static int sy_get_parenthesis(char chr)
{
    if (chr == '(')
        return PARENTHESIS_LEFT;
    if (chr == ')')
        return PARENTHESIS_RIGHT;
    return PARENTHESIS_NONE;
}

/**
 * Helper function to retrieve function token identifier
 * if not a function, return -1 (FUNC_UNSUPPORTED constant)
 */
static int sy_get_function(char** chr)
{
    int i;
    int size = (int) (sizeof(func_match) / sizeof(func_match_template));

    /* iterate through supported function map and try to match at least one of them */
    for (i = 0; i < size; i++)
    {
        /* compares function name with string on input */
        if (strncmp(*chr, func_match[i].func_name, strlen(func_match[i].func_name)) == 0)
        {
            /* if matched, move the input pointer by function identifier length */
            *chr += strlen(func_match[i].func_name);
            return func_match[i].func_id;
        }
    }

    return FUNC_UNSUPPORTED;
}

/**
 * Helper function to retrieve variable identifier
 * this method is highly customized to specification of semestral work
 * - parses only one character string and returns positive value as "found" if the char
 *   is alphanumeric
 */
static int sy_get_variable(char** chr)
{
    char var = **chr;

    /* parse only one-character (one letter) variables */
    if ((var >= 'a' && var <= 'z') || (var >= 'A' && var <= 'Z'))
    {
        *chr += 1;
        return (int)var;
    }

    return -1;
}

/*
 * Compares two operators and decides the priority (has greater, equal or lower)
 * returns 1 if first has greater priority than second (i.e. * over +)
 * returns 0 if both operators are equal
 * returns -1 if first has lower priority than second (i.e. + and *)
 */
static int sy_operator_priority(int first, int second)
{
    /* exponential operator ^ has always greater priority than other operators */

    if (first == OP_EXP_RAISE && second == OP_EXP_RAISE)
        return 0;

    if (first == OP_EXP_RAISE)
        return 1;

    if (second == OP_EXP_RAISE)
        return -1;

    /* test different "priority classes" */

    if ((first == OP_ADD || first == OP_SUBTRACT) && (second == OP_MULTIPLY || second == OP_DIVIDE))
        return -1;

    if ((first == OP_MULTIPLY || first == OP_DIVIDE) && (second == OP_ADD || second == OP_SUBTRACT))
        return 1;

    /* test same "priority classes" */

    if ((first == OP_ADD || first == OP_SUBTRACT) && (second == OP_ADD || second == OP_SUBTRACT))
        return 0;

    if ((first == OP_MULTIPLY || first == OP_DIVIDE) && (second == OP_MULTIPLY || second == OP_DIVIDE))
        return 0;

    /* now everything should be tested, so anything else is considered error */

    return 2;
}

/**
 * Generates RPN represented expression from input string
 * if something fails, sets flag and returns position of character, where everything failed
 */
c_stack* sy_generate_rpn_stack(char *input, int *error, char** error_ptr)
{
    /* output stack of rpn_elements */
    c_stack *rpn_stack;
    /* temporary stack of operators */
    c_stack *op_stack;

    rpn_element *rpn_el_tmp;

    char chr;
    int tmp, sign, stored_sign, flag, digit;
    double dtmp, weight;
    rpn_element *last_el;

    if (strlen(input) == 0)
    {
        *error = SYNTAX_ERROR_NOTHING_TO_PARSE;
        return NULL;
    }

    /* create stacks used in Shunting yard algorhitm */
    rpn_stack = stck_create(RPN_STACK_SIZE);
    op_stack = stck_create(SY_OP_STACK_SIZE);

    *error = SYNTAX_ERROR_NONE;
    *error_ptr = NULL;
    stored_sign = 1;
    last_el = NULL;

    /* go through string character by character */
    while (*input != (char)0)
    {
        chr = *input;

        /* first, try to parse numeric constant (integer / real) */
        tmp = sy_get_number(chr);
        /* note that "dot" character is also considered as beginning of numeric constant */
        if (tmp != -1 || chr == '.')
        {
            /* if the function is followed by constant, not parenthesis, it's error */
            if (last_el != NULL && last_el->type == RPN_TOKEN_FUNCTION)
            {
                ERROR_ROUTINE(SYNTAX_ERROR_FUNCTION_PARENTHESIS, input);
                return NULL;
            }

            /* this allows parsing of real numbers in i.e. ".5E2" format */
            /* if the dot is not the first character... */
            if (chr != '.')
            {
                /* parse the integer part as it should be */
                while ((digit = sy_get_number(*(++input))) != (int)-1)
                {
                    tmp = tmp * 10 + digit;
                }
            }
            /* if it begins with dot, the integer part is equal 0 */
            else
            {
                tmp = 0;
            }

            dtmp = (double)(tmp);

            /* if non-integer part is present */
            if (*input == '.')
            {
                /* parse it and add to working value */
                weight = 0.1;
                while ((digit = sy_get_number(*(++input))) != (int)-1)
                {
                    dtmp += (double)(digit)*weight;
                    weight /= 10.0;
                }

                /* if decadic exponent is present */
                if (*input == 'E' || *input == 'e')
                {
                    sign = 1;
                    tmp = 0;

                    /* move pointer (so we won't parse 'e' again) */
                    chr = *(++input);

                    /* the exponent may have sign */
                    if (chr == '-')
                    {
                        sign = -1;
                    }
                    else if (chr == '+')
                    {
                        /* sign already set to 1 */
                    }
                    else
                    {
                        /* the exponent mark ('e' or 'E') has to be followed by number */
                        if ((digit = sy_get_number(chr)) != (int)-1)
                            tmp = digit;
                        else
                        {
                            ERROR_ROUTINE(SYNTAX_ERROR_REAL_NOTATION, input);
                            return NULL;
                        }
                    }

                    /* parse the value of exponent */
                    while ((digit = sy_get_number(*(++input))) != (int)-1)
                    {
                        tmp  = tmp*10 + digit;
                    }

                    /* and proceed it to working value */
                    dtmp = dtmp * pow(10.0, (double)(sign*tmp));
                }
            }

            /* apply stored sign (may be stored i.e. from parsing exponent operator */
            dtmp = dtmp*stored_sign;
            stored_sign = 1;

            /* build element and push it to stack */
            rpn_el_tmp = rpn_build_element(RPN_TOKEN_CONST);
            last_el = rpn_el_tmp;
            rpn_el_tmp->value.as_double = dtmp;
            stck_push(rpn_stack, rpn_el_tmp);

            continue;
        }

        /* try to parse operator */
        tmp = sy_get_operator(chr);
        if (tmp != OP_NONE)
        {
            /* if the function is followed by operator, not parenthesis, it's error */
            if (last_el != NULL && last_el->type == RPN_TOKEN_FUNCTION)
            {
                ERROR_ROUTINE(SYNTAX_ERROR_FUNCTION_PARENTHESIS, input);
                return NULL;
            }

            /* let's validate unary operators and operands of binary operators */
            if ((last_el == NULL && rpn_stack->curr == STCK_INVALID) || (last_el != NULL && last_el->type == RPN_TOKEN_OPERATOR && last_el->value.as_operator != PARENTHESIS_RIGHT))
            {
                /* unary operators + and - */
                if ((tmp == OP_ADD || tmp == OP_SUBTRACT))
                {
                    /* if two same-priority operators came, that's pretty suspicious */
                    if (last_el != NULL && (last_el->value.as_operator == OP_ADD || last_el->value.as_operator == OP_SUBTRACT))
                    {
                        /* unary minus can appear directly after binary plus - anything else in this block is considered error */
                        if (tmp != OP_SUBTRACT || last_el->value.as_operator != OP_ADD)
                        {
                            ERROR_ROUTINE(SYNTAX_ERROR_OPERATOR_FREQUENCY, input);
                            return NULL;
                        }
                    }

                    /* unary operator after exponent mark */
                    if (last_el != NULL && last_el->value.as_operator == OP_EXP_RAISE)
                    {
                        if (tmp == OP_SUBTRACT)
                        {
                            stored_sign = -1;
                        }

                        input++;
                        continue;
                    }

                    /* insert "hidden zero" element, because -5 (unary minus) has the same value as 0-5, as well, as +5 and 0+5 */
                    rpn_el_tmp = rpn_build_element(RPN_TOKEN_CONST);
                    last_el = rpn_el_tmp;
                    rpn_el_tmp->value.as_double = 0;
                    stck_push(rpn_stack, rpn_el_tmp);
                }
                else
                {
                    ERROR_ROUTINE(SYNTAX_ERROR_BINARY_OPERATOR_OPERANDS, input);
                    return NULL;
                }
            }

            /* "look" at last element in operator stack */
            rpn_el_tmp = (rpn_element*)stck_peek(op_stack);

            if (rpn_el_tmp != NULL)
            {
                sign = sy_operator_priority(tmp, rpn_el_tmp->value.as_operator);
                /* lower priority operator came */
                if (sign == -1 || (sign == 0 && tmp != OP_EXP_RAISE))
                {
                    /* "eat" all greater or equal priority operators and push them onto RPN stack */
                    while ((rpn_el_tmp = (rpn_element*)stck_peek(op_stack)) != NULL)
                    {
                        if (sy_operator_priority(tmp, rpn_el_tmp->value.as_operator) <= 0)
                        {
                            stck_push(rpn_stack, stck_pop(op_stack));
                        }
                        else
                        {
                            break;
                        }
                    }
                }
                /* greater priority operator came */
                else if (sign == 1)
                {
                    /* just push operator element on stack - handled below */
                }
                /* exponential operator came, and previous also was exponential (special case, only one right-associative operator) */
                else if (sign == 0 && tmp == OP_EXP_RAISE)
                {
                    /* just push operator element on stack - handled below */
                }
            }

            rpn_el_tmp = rpn_build_element(RPN_TOKEN_OPERATOR);
            last_el = rpn_el_tmp;
            rpn_el_tmp->value.as_operator = tmp;
            stck_push(op_stack, rpn_el_tmp);
            ++input;
            continue;
        }

        /* now try to parse parenthesis */
        tmp = sy_get_parenthesis(chr);
        if (tmp != PARENTHESIS_NONE)
        {
            /* left parenthesis is simply pushed to operator stack */
            if (tmp == PARENTHESIS_LEFT)
            {
                rpn_el_tmp = rpn_build_element(RPN_TOKEN_OPERATOR);
                last_el = rpn_el_tmp;
                rpn_el_tmp->value.as_operator = tmp;
                stck_push(op_stack, rpn_el_tmp);
            }
            else
            {
                /* closing bracket immediatelly after function token is also wrong */
                if (last_el != NULL && last_el->type == RPN_TOKEN_FUNCTION)
                {
                    ERROR_ROUTINE(SYNTAX_ERROR_FUNCTION_PARENTHESIS, input);
                    return NULL;
                }

                flag = -1;

                /* go through operators in operator stack and push operators, that aren't parenthesis
                   on rpn stack, until we hit left parenthesis */
                while ((rpn_el_tmp = (rpn_element*)stck_pop(op_stack)) != NULL)
                {
                    if (rpn_el_tmp->value.as_operator == PARENTHESIS_LEFT)
                    {
                        /* parenthesis found flag */
                        flag = 1;

                        if (last_el == rpn_el_tmp)
                            last_el = NULL;

                        /* we won't reuse this, let's free it */
                        free(rpn_el_tmp);

                        /* also look, if parentheses was used to match math function argument */
                        rpn_el_tmp = (rpn_element*)stck_peek(op_stack);
                        if (rpn_el_tmp != NULL && rpn_el_tmp->type == RPN_TOKEN_FUNCTION)
                            stck_push(rpn_stack, stck_pop(op_stack));

                        break;
                    }

                    last_el = NULL;
                    stck_push(rpn_stack, rpn_el_tmp);
                }

                /* if we hit the bottom without parenthesis being found, it means that opening parenthesis is missing */
                if (flag == -1)
                {
                    ERROR_ROUTINE(SYNTAX_ERROR_MISSING_PARENTHESIS, input);
                    return NULL;
                }
            }

            ++input;
            continue;
        }

        /* now try to parse function */
        tmp = sy_get_function(&input);
        if (tmp != FUNC_UNSUPPORTED)
        {
            /* functions are just thrown onto stack */
            rpn_el_tmp = rpn_build_element(RPN_TOKEN_FUNCTION);
            last_el = rpn_el_tmp;
            rpn_el_tmp->value.as_function = tmp;
            stck_push(op_stack, rpn_el_tmp);
            continue;
        }

        /* as last thing, try to parse variable (one letter token) */
        tmp = sy_get_variable(&input);
        if (tmp != -1)
        {
            /* as well, as constant, the variable should not appear as next element after function
             * token without parenthesis being immediatelly open */
            if (last_el != NULL && last_el->type == RPN_TOKEN_FUNCTION)
            {
                ERROR_ROUTINE(SYNTAX_ERROR_FUNCTION_PARENTHESIS, input);
                return NULL;
            }

            rpn_el_tmp = rpn_build_element(RPN_TOKEN_VARIABLE);
            last_el = rpn_el_tmp;
            rpn_el_tmp->value.as_variable = tmp;
            stck_push(rpn_stack, rpn_el_tmp);
            continue;
        }

        /* ignore whitespaces */
        if (chr == ' ')
        {
            ++input;
            continue;
        }

        /* this means we have some unexpected character / token on input */
        ERROR_ROUTINE(SYNTAX_ERROR_INVALID_CHARACTER, input);
        return NULL;
    }

    /* if last element of expression was an operator (binary), it's an error */
    if (last_el != NULL && last_el->type == RPN_TOKEN_OPERATOR)
    {
        ERROR_ROUTINE(SYNTAX_ERROR_BINARY_OPERATOR_OPERANDS, input);
        return NULL;
    }

    /* every operator left in operator stack should be pushed to output stack */
    while ((rpn_el_tmp = (rpn_element*)stck_pop(op_stack)) != NULL)
    {
        /* ..except parenthesis
         * this also validates proper closing of parenthesis */
        if (rpn_el_tmp->type == RPN_TOKEN_OPERATOR && rpn_el_tmp->value.as_operator == PARENTHESIS_LEFT)
        {
            free(rpn_el_tmp);
            ERROR_ROUTINE(SYNTAX_ERROR_MISSING_PARENTHESIS, input);
            return NULL;
        }
        stck_push(rpn_stack, rpn_el_tmp);
    }

    /* clean up temporary operator stack */
    stck_destroy(op_stack);

    return rpn_stack;
}
