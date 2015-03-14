#ifndef MATHPARSER_STACK_H
#define MATHPARSER_STACK_H

#define STCK_INVALID -1

/* "c" stands for "custom" */
struct c_stack
{
    void **elements;
    int curr;
    int size;
};

c_stack* stck_create(int size);
void stck_destroy(c_stack *stck);
void* stck_pop(c_stack *stck);
void* stck_peek(c_stack *stck);
void* stck_get(c_stack *stck, int pos);
void stck_push(c_stack *stck, void *el);

#endif
