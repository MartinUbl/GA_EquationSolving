#ifndef MATHPARSER_RPN_H
#define MATHPARSER_RPN_H

enum rpn_token_type
{
    RPN_TOKEN_CONST,                /* all numeric constants */
    RPN_TOKEN_VARIABLE,             /* variable(s) */
    RPN_TOKEN_OPERATOR,             /* binary operators */
    RPN_TOKEN_FUNCTION              /* math functions */
};

typedef struct
{
    enum rpn_token_type type;
    union
    {
        int as_variable;                        /* for variables */
        double as_double;                       /* for constants */
        int as_operator;                        /* for operators */
        int as_function;                        /* for recognized math functions */
    } value;
} rpn_element;

rpn_element* rpn_build_element(enum rpn_token_type type);
double rpn_evaluate_stack(c_stack* stck, double* variable_values);

#endif
