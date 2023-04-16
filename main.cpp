#include <md_dry.h>

int main(int argc, char const *argv[]){
    MDdry *md;
    md = new MDdry();
    md->make_trajectory();
    return 0;
}