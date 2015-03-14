#ifndef MATHPARSER_MAIN_H
#define MATHPARSER_MAIN_H

/*
 * MACROS
 */

#define RPN_STACK_SIZE 64
#define SY_OP_STACK_SIZE 32

#define PLOT_STEP_COEF 0.005
#define PLOT_SIZE_X 200
#define PLOT_SIZE_Y 200

#ifndef M_PI /* i.e. MSVS case */
#define M_PI 3.14159265
#endif

/*
 * ENUMS
 */

enum operator_type
{
    OP_NONE = -1,                   /* error flag */
    OP_ADD,                         /*   +   */
    OP_SUBTRACT,                    /*   -   */
    OP_MULTIPLY,                    /*   *   */
    OP_DIVIDE,                      /*   /   */
    OP_EXP_RAISE,                   /*   ^   */
    /* not really operators, but sy parser deals with them like that */
    PARENTHESIS_LEFT,               /*   (   */
    PARENTHESIS_RIGHT               /*   )   */
};

#define PARENTHESIS_NONE OP_NONE    /* error flag to allow generic handling with disambiguation for programmer */

enum supported_functions
{
    FUNC_UNSUPPORTED = -1,          /* error flag */

    FUNC_ABS,                       /* absolute value */
    FUNC_EXP,                       /* exponential, e^x */

    FUNC_SIN,                       /* sinus */
    FUNC_COS,                       /* cosinus */
    FUNC_TAN,                       /* tangens */
    FUNC_COTAN,                     /* cotangens */

    FUNC_ASIN,                      /* arcus sinus */
    FUNC_ACOS,                      /* arcus cosinus */
    FUNC_ATAN,                      /* arcus tangens */
    FUNC_ACOTAN,                    /* arcus cotangens */

    FUNC_LOG10,                     /* logarhitm base 10 */
    FUNC_LN,                        /* logarhitm base e */

    FUNC_SINH,                      /* hyperbolic sinus */
    FUNC_COSH,                      /* hyperbolic cosinus */
    FUNC_TANH,                      /* hyperbolic tangens */

    FUNC_TODEG,                     /* convert radians to degrees */
    FUNC_TORAD                      /* convert degrees to radians */
};

enum syntax_error_code
{
    SYNTAX_ERROR_NONE = 0,
    SYNTAX_ERROR_MISSING_PARENTHESIS,
    SYNTAX_ERROR_REAL_NOTATION,
    SYNTAX_ERROR_OPERATOR_FREQUENCY,
    SYNTAX_ERROR_FUNCTION_PARENTHESIS,
    SYNTAX_ERROR_BINARY_OPERATOR_OPERANDS,
    SYNTAX_ERROR_INVALID_CHARACTER,
    SYNTAX_ERROR_NOTHING_TO_PARSE
};

#endif
