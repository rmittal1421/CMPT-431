#include "../common/allocator.h"
#include "../common/utils.h"

#define LFENCE asm volatile("lfence" : : : "memory")
#define SFENCE asm volatile("sfence" : : : "memory")

template<class P>
class pointer_t {
    public:
        P* ptr;

        P* address() {
            ptr = (Node<T>*) ((uintptr_t)ptr << 16);
            ptr = (Node<T>*) (uintptr_t)ptr >> 16;
            return ptr;
        }

        uint count() {
            uint count = (uint) (uintptr_t)ptr >> 48;
            return count;
        }
};

template <class T>
class Node
{
public:
    T value;
    pointer_t<Node<T>> next;
};

template <class T>
class NonBlockingQueue
{
    pointer_t<Node<T>> q_head;
    pointer_t<Node<T>> q_tail;
    CustomAllocator my_allocator_;
public:
    
    NonBlockingQueue() : my_allocator_()
    {
        std::cout << "Using NonBlockingQueue\n";
    }

    void initQueue(long t_my_allocator_size){
        std::cout << "Using Allocator\n";
        my_allocator_.initialize(t_my_allocator_size, sizeof(Node<T>));
        // Initialize the queue head or tail here
        Node<T>* newNode = (Node<T>*)my_allocator_.newNode();
        newNode->next.ptr = nullptr;
        q_head.ptr = newNode;
        q_tail.ptr = newNode;
        // my_allocator_.freeNode(newNode);
    }

    void enqueue(T value)
    {
        // Use LFENCE and RFENCE as mentioned in pseudocode

        Node<T>* node = (Node<T>* )my_allocator_.newNode();
        node->value = value;
        node->next.ptr = nullptr;

        SFENCE;
        pointer_t<Node<T>> tail;

        while(true) {
            tail = q_tail;
            LFENCE;
            pointer_t<Node<T>> next = tail.address()->next;
            LFENCE;
            if(tail.ptr == q_tail.ptr) {
                if(next.address() == nullptr) {
                    // CAS operation
                    if(CAS(&tail.address()->next, next, <node, next.count() + 1>)) {
                        break;
                    }
                } else {
                    // CAS operation
                    CAS(&q_tail, tail, <next.address(), tail.count() + 1>);
                }
            }
        }

        SFENCE;
        // CAS operation
        CAS(&q_tail, tail, <node, tail.count() + 1>);
    }

    bool dequeue(T *value)
    {
        // Use LFENCE and RFENCE as mentioned in pseudocode

        pointer_t<Node<T>> head;
        
        while(true) {
            head = q_head;
            LFENCE;
            pointer_t<Node<T>> tail = q_tail;
            LFENCE;
            pointer_t<Node<T>> next = head.address()->next;
            LFENCE;

            if(head.ptr == q_head.ptr) {
                if(head.address() == tail.address()) {
                    if(next.address() == nullptr) {
                        return false;
                    }
                    CAS(&q_tail, tail, <next.address(), tail.count() + 1>);
                } else {
                    *value = next.address()->value;
                    if(CAS(&q_head, head, <next.address(), head.count() + 1>)) {
                        break;
                    }
                }
            }
        }

        my_allocator_.freeNode(head.address());
        return true;
    }

    void cleanup()
    {
        my_allocator_.cleanup();
    }

};

