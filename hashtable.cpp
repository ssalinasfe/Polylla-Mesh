#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "hashtable.h"


#ifdef DEBUG
#define DEBUG_TEST 1
#else
#define DEBUG_TEST 0
#endif

#define debug_block(fmt) do { if (DEBUG_TEST){ fmt }} while (0)
#define debug_print(fmt, ...) do { if (DEBUG_TEST) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, __LINE__, __func__, __VA_ARGS__); } while (0)
#define debug_msg(fmt) do { if (DEBUG_TEST) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__,  __LINE__, __func__); } while (0)

int MAX_HASH = 200000;



int hashy(int key){
	return key % MAX_HASH;
}

void generate_hash_table(node** hashtable){
    for(int i = 0; i<MAX_HASH; i++){
        hashtable[i] = NULL;
    }
}

//crea un nodo
node* create(int triangle, node* next){
    node* new_node = (node*)malloc(sizeof(node));
    if(new_node == NULL){
        printf("Error creating a new node.\n");
        exit(0);
    }
    new_node->triangle = triangle;
    new_node->next = next;
    return new_node;
}

// agrega un nodo al inicio de la lista
node* prepend(node* head, int triangle){
    node* new_node = create(triangle, head);
    head = new_node;
    return head;
}



//muestra los elementos de la lista dinamica
void mostrar(node* head){
    node* cursor = head;
    while(cursor != NULL){
        std::cout<<cursor->triangle<<std::endl;
        cursor = cursor->next;
    }
}

//busca un elemento en el nodo
int search(node* head,int triangle){
    node *cursor = head;
    while(cursor!=NULL){
        if(cursor->triangle == triangle)
            return 1;
        cursor = cursor->next;
    }
    return 0;
}


// remove an element from the linked list
void searchandremove(node *&head, int triangle){
    node* current = head; // the first valid node
    node* prev = NULL; // empty header
    debug_print("Revisando triangulo %d en hashtable\n", triangle);
    if(current != NULL && current->triangle == triangle){
        debug_print("Triangulo %d repetido en HEAD\n", current->triangle);
        head = current->next;
        free(current);
    }else{
        while(current != NULL && current->triangle != triangle){
            debug_msg("dasd");
            prev = current; 
            current = current->next; // go to next value
        }
        if(current!=NULL){
            debug_print("Triangulo %d repetido, borrando por %d\n", triangle, current->triangle);
            prev->next = current->next;
            free(current);
        }
    } 
}



//busca un elemento en el nodo
//ZZ return_distancia(node* head,ZZ trampa){
//    node *cursor = head;
//    while(cursor!=NULL){
//        if(cursor->trampa == trampa)
//            return cursor->distancia;
//        cursor = cursor->next;
//    }
//    return conv<ZZ>(0);
// }

//borra la lista
void dispose(node *head){
    node *cursor, *tmp;
    if(head != NULL){
        cursor = head->next;
        head->next = NULL;
        while(cursor != NULL){
            tmp = cursor->next;
            free(cursor);
            cursor = tmp;
        }
    }
}