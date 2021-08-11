extern int MAX_HASH;

typedef struct node{
    int triangle;
    struct node* next;
} node;

node* prepend(node* head, int triangle); // insterta un elemento al inicio de la lista
node* create(int triangle, node* next); //crea un nodo nuevo
void mostrar(node* head); // muestra todos los elementos de la lista
int search(node* head,int triangle); // busca un elemento en la lista
void dispose(node *head); // borra la lista
//ZZ return_distancia(node* head,ZZ trampa);
int hashy(int key);
void generate_hash_table(node** hashtable);
void searchandremove(node *&head, int triangle);