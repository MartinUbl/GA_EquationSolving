#ifndef MATHPARSER_SY_H
#define MATHPARSER_SY_H

c_stack* sy_generate_rpn_stack(char *input, int *error, char** error_ptr);

struct _func_match_template
{
    int func_id;
    const char* func_name;
};

typedef struct _func_match_template func_match_template;

#endif
