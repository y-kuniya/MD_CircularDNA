#include <md.h>

int main(int argc, char const *argv[]){
    MD *md;
    md = new MD();
    md->make_trajectory();
    return 0;
}