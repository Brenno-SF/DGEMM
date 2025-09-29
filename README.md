# DGEMM Sequencial e Paralelo com OpenMP

Este projeto implementa e compara a multiplicação de matrizes (`DGEMM`) de forma **sequencial** e **paralela** utilizando **OpenMP**.  
São aplicadas técnicas de otimização como **tiling (blocking)**, **transposição da matriz B** e **register blocking** para melhor aproveitamento de cache.

---


**Universidade:** Universidade Estadual de Santa Cruz  
**Curso:** Engenharia de Computação / Processamento Paralelo  
**Disciplina:** DEC107 - Processamento Paralelo  
**Professor:** Dr. Esbel Tomás Valero Orellana  
**Alunos:** Brenno S. Florêncio e Mateus Soares  
**Data:** 28/09/2025 

---

## 🔧 Compilação

Compile o código no Linux com o comando:

```bash
gcc -o dgemmSeqPar -Wall -O3 -fopenmp -march=native -mfma dgemmSeqParGPT.c -lm
```

Explicação das flags:  
- `-Wall` → mostra todos os avisos do compilador.  
- `-O3` → nível alto de otimização.  
- `-fopenmp` → habilita OpenMP.  
- `-march=native` → usa instruções específicas da sua CPU.  
- `-mfma` → habilita instruções FMA (fused multiply-add) se disponíveis.  
- `-lm` → linka a biblioteca matemática.

---

## ▶️ Execução

Rode o programa passando o tamanho da matriz quadrada `N`:

```bash
./dgemmSeqPar
Matrix size: 1024
```

---

## 📊 Saída esperada

O programa gera como saída:

1. **Tempo de execução** da versão sequencial e das versões paralelas com 2, 4, 8 e 12 threads.  
2. **Speedup** → quanto a versão paralela foi mais rápida em relação à sequencial.  
3. **Eficiência** → quão bem os threads foram aproveitados (`speedup / número de threads`).  
4. **Diferença numérica** entre os resultados sequencial e paralelo (verificação de corretude).

Exemplo de saída:

```
dgemmSeq time: 0.337004 seconds
dgemmPar with 2 threads time: 0.188320 seconds
dgemmPar with 4 threads time: 0.126693 seconds
dgemmPar with 8 threads time: 0.122472 seconds
dgemmPar with 12 threads time: 0.119871 seconds

Speedup with 2 threads: 1.78
Speedup with 4 threads: 2.66
Speedup with 8 threads: 2.75
Speedup with 12 threads: 2.81

Efficiency with 2 threads: 0.89
Efficiency with 4 threads: 0.66
Efficiency with 8 threads: 0.34
Efficiency with 12 threads: 0.23

Calculating differences between sequential and parallel results...
Difference with 2 threads: 0.000000e+00
Difference with 4 threads: 0.000000e+00
Difference with 8 threads: 0.000000e+00
Difference with 12 threads: 0.000000e+00
```

---

## 📚 Técnicas usadas

- **Transposição de B** → melhora a localidade de acesso.  
- **Blocking (tiling)** → divide as matrizes em blocos de 64x64 para melhor uso do cache.  
- **Register blocking** → acumula parciais em registradores para reduzir acessos à memória.  
- **OpenMP** → paralelização do loop externo com `#pragma omp parallel for collapse(2)`.

---
