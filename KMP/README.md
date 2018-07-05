compile: gcc gen_ped.c ¨Co gen_ped
¡¡¡¡¡¡
         mpicc kmp.c ¨Co kmp


execute£º
        First, execute gen_ped to generate pattern string,
        *gen_ped Strlen Pedlen Seed Pattern_File*

        **Strlen** is the length of the pattern string£¬
        **Pedlen** is the minimum period length of the pattern string
        **Seed** is used for random functions
        **Pattern_File** is the file that includes the generated data, here in 
        kmp.c it is fixed as "pattern.dat"

        Then, execute kmp