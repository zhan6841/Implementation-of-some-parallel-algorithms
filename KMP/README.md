compile:   
&emsp;&emsp;&emsp;  gcc gen_ped.c -o gen_ped   
&emsp;&emsp;&emsp;  mpicc kmp.c -o kmp          
execute:   
&emsp;&emsp;&emsp;  First, execute gen_ped to generate pattern string   
&emsp;&emsp;&emsp;&emsp;  *gen_ped Strlen Pedlen Seed Pattern_File*   
&emsp;&emsp;&emsp;&emsp;  **Strlen** is the length of the pattern string   
&emsp;&emsp;&emsp;&emsp;  **Pedlen** is the minimum period length of the pattern string   
&emsp;&emsp;&emsp;&emsp;  **Seed** is used for random functions   
&emsp;&emsp;&emsp;&emsp;  **Pattern_File** is the file that includes the generated data, here in kmp.c it is fixed as "pattern.dat"     
&emsp;&emsp;&emsp;  Then, execute kmp
