#include "../common/allocator.h"
#include <mutex>

template <class T>
class Node
{
    T value;
    Node<T>* next;
};

template <class T>
class OneLockQueue
{
    Node<T>* q_head;
    Node<T>* q_tail;
    CustomAllocator my_allocator_;
    std::mutex lock;

public:
    OneLockQueue() : my_allocator_()
    {
        std::cout << "Using OneLockQueue\n";
    }

    void initQueue(long t_my_allocator_size){
        std::cout << "Using Allocator\n";
        my_allocator_.initialize(t_my_allocator_size, sizeof(Node<T>));
        // Initialize the queue head or tail here
        Node<T>* newNode = (Node<T>*)my_allocator_.newNode();
        newNode->next = nullptr;
        q_head = newNode;
        q_tail = newNode;
        // my_allocator_.freeNode(newNode);
    }

    void enqueue(T value)
    {
        Node<T>* newNode = (Node<T>*)my_allocator_.newNode();
        node->value = value;
        node->next = nullptr;

        lock.lock();
        q_tail->next = node;
        q_tail = node;
        lock.unlock();
    }

    bool dequeue(T *value)
    {
        lock.lock();
        Node<T>* node = q_head;
        Node<T>* new_head = q_head->next;
        if(new_head == nullptr) {
            lock.unlock();
            return false;
        }
        *value = new_head->value;
        q_head = new_head;
        my_allocator_.freeNode(node);
        lock.unlock();
        return ret_value;
    }

    void cleanup()
    {
        my_allocator_.cleanup();
    }
};