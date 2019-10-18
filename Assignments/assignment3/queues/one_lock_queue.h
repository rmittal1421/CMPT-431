#include "../common/allocator.h"

template <class T>
class Node
{
};

template <class T>
class OneLockQueue
{
    CustomAllocator my_allocator_;
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
        my_allocator_.freeNode(newNode);
    }

    void enqueue(T value)
    {
    }

    bool dequeue(T *value)
    {
        bool ret_value = false;
        return ret_value;
    }

    void cleanup()
    {
        my_allocator_.cleanup();
    }
};