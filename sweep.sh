#!/usr/bin/env bash
set -e

# Make outdirput folders
mkdir -p runs/served runs/time

# Make outdirput folders
mkdir -p runs/served runs/time

# ========== BY_SERVED (N = 100,000 customers) ==========
# Light/medium/heavy traffic mixes; vary μ at fixed λ, and vary λ at fixed μ.

# Fixed λ = 0.70; vary μ
./des_sim --lambda 0.70 --mu 0.90 --term served --max-served 100000 --seed 11 --warmup 1000 --reps 10 --outdir runs/served/l070_m090_s11/
./des_sim --lambda 0.70 --mu 1.00 --term served --max-served 100000 --seed 11 --warmup 1000 --reps 10 --outdir runs/served/l070_m100_s11/
./des_sim --lambda 0.70 --mu 1.10 --term served --max-served 100000 --seed 11 --warmup 1000 --reps 10 --outdir runs/served/l070_m110_s11/
./des_sim --lambda 0.70 --mu 1.20 --term served --max-served 100000 --seed 11 --warmup 1000 --reps 10 --outdir runs/served/l070_m120_s11/
./des_sim --lambda 0.70 --mu 1.30 --term served --max-served 100000 --seed 11 --warmup 1000 --reps 10 --outdir runs/served/l070_m130_s11/
./des_sim --lambda 0.70 --mu 1.40 --term served --max-served 100000 --seed 11 --warmup 1000 --reps 10 --outdir runs/served/l070_m140_s11/
./des_sim --lambda 0.70 --mu 1.50 --term served --max-served 100000 --seed 11 --warmup 1000 --reps 10 --outdir runs/served/l070_m150_s11/
./des_sim --lambda 0.70 --mu 1.60 --term served --max-served 100000 --seed 11 --warmup 1000 --reps 10 --outdir runs/served/l070_m160_s11/
./des_sim --lambda 0.70 --mu 1.70 --term served --max-served 100000 --seed 11 --warmup 1000 --reps 10 --outdir runs/served/l070_m170_s11/
./des_sim --lambda 0.70 --mu 1.80 --term served --max-served 100000 --seed 11 --warmup 1000 --reps 10 --outdir runs/served/l070_m180_s11/

# Fixed μ = 1.10; vary λ
./des_sim --lambda 0.10 --mu 1.10 --term served --max-served 100000 --seed 11 --warmup 1000 --reps 10 --outdir runs/served/l010_m110_s11/
./des_sim --lambda 0.20 --mu 1.10 --term served --max-served 100000 --seed 11 --warmup 1000 --reps 10 --outdir runs/served/l020_m110_s11/
./des_sim --lambda 0.30 --mu 1.10 --term served --max-served 100000 --seed 11 --warmup 1000 --reps 10 --outdir runs/served/l030_m110_s11/
./des_sim --lambda 0.40 --mu 1.10 --term served --max-served 100000 --seed 11 --warmup 1000 --reps 10 --outdir runs/served/l040_m110_s11/
./des_sim --lambda 0.50 --mu 1.10 --term served --max-served 100000 --seed 11 --warmup 1000 --reps 10 --outdir runs/served/l050_m110_s11/
./des_sim --lambda 0.60 --mu 1.10 --term served --max-served 100000 --seed 11 --warmup 1000 --reps 10 --outdir runs/served/l060_m110_s11/
./des_sim --lambda 0.70 --mu 1.10 --term served --max-served 100000 --seed 11 --warmup 1000 --reps 10 --outdir runs/served/l070_m110_s11/
./des_sim --lambda 0.80 --mu 1.10 --term served --max-served 100000 --seed 11 --warmup 1000 --reps 10 --outdir runs/served/l080_m110_s11/
./des_sim --lambda 0.90 --mu 1.10 --term served --max-served 100000 --seed 11 --warmup 1000 --reps 10 --outdir runs/served/l090_m110_s11/
./des_sim --lambda 1.00 --mu 1.10 --term served --max-served 100000 --seed 11 --warmup 1000 --reps 10 --outdir runs/served/l100_m110_s11/

# ========== BY_TIME (T = 50,000 time units) ==========
# Light/medium/heavy traffic mixes; vary μ at fixed λ, and vary λ at fixed μ.


# Fixed λ = 0.70; vary μ
./des_sim --lambda 0.70 --mu 0.90 --term time --horizon 50000 --seed 11 --warmup 1000 --reps 10 --outdir runs/time/l070_m090_s11/
./des_sim --lambda 0.70 --mu 1.00 --term time --horizon 50000 --seed 11 --warmup 1000 --reps 10 --outdir runs/time/l070_m100_s11/
./des_sim --lambda 0.70 --mu 1.10 --term time --horizon 50000 --seed 11 --warmup 1000 --reps 10 --outdir runs/time/l070_m110_s11/
./des_sim --lambda 0.70 --mu 1.20 --term time --horizon 50000 --seed 11 --warmup 1000 --reps 10 --outdir runs/time/l070_m120_s11/
./des_sim --lambda 0.70 --mu 1.30 --term time --horizon 50000 --seed 11 --warmup 1000 --reps 10 --outdir runs/time/l070_m130_s11/
./des_sim --lambda 0.70 --mu 1.40 --term time --horizon 50000 --seed 11 --warmup 1000 --reps 10 --outdir runs/time/l070_m140_s11/
./des_sim --lambda 0.70 --mu 1.50 --term time --horizon 50000 --seed 11 --warmup 1000 --reps 10 --outdir runs/time/l070_m150_s11/
./des_sim --lambda 0.70 --mu 1.60 --term time --horizon 50000 --seed 11 --warmup 1000 --reps 10 --outdir runs/time/l070_m160_s11/
./des_sim --lambda 0.70 --mu 1.70 --term time --horizon 50000 --seed 11 --warmup 1000 --reps 10 --outdir runs/time/l070_m170_s11/
./des_sim --lambda 0.70 --mu 1.80 --term time --horizon 50000 --seed 11 --warmup 1000 --reps 10 --outdir runs/time/l070_m180_s11/

# Fixed μ = 1.10; vary λ
./des_sim --lambda 0.10 --mu 1.10 --term time --horizon 50000 --seed 11 --warmup 1000 --reps 10 --outdir runs/time/l010_m110_s11/
./des_sim --lambda 0.20 --mu 1.10 --term time --horizon 50000 --seed 11 --warmup 1000 --reps 10 --outdir runs/time/l020_m110_s11/
./des_sim --lambda 0.30 --mu 1.10 --term time --horizon 50000 --seed 11 --warmup 1000 --reps 10 --outdir runs/time/l030_m110_s11/
./des_sim --lambda 0.40 --mu 1.10 --term time --horizon 50000 --seed 11 --warmup 1000 --reps 10 --outdir runs/time/l040_m110_s11/
./des_sim --lambda 0.50 --mu 1.10 --term time --horizon 50000 --seed 11 --warmup 1000 --reps 10 --outdir runs/time/l050_m110_s11/
./des_sim --lambda 0.60 --mu 1.10 --term time --horizon 50000 --seed 11 --warmup 1000 --reps 10 --outdir runs/time/l060_m110_s11/
./des_sim --lambda 0.70 --mu 1.10 --term time --horizon 50000 --seed 11 --warmup 1000 --reps 10 --outdir runs/time/l070_m110_s11/
./des_sim --lambda 0.80 --mu 1.10 --term time --horizon 50000 --seed 11 --warmup 1000 --reps 10 --outdir runs/time/l080_m110_s11/
./des_sim --lambda 0.90 --mu 1.10 --term time --horizon 50000 --seed 11 --warmup 1000 --reps 10 --outdir runs/time/l090_m110_s11/
./des_sim --lambda 1.00 --mu 1.10 --term time --horizon 50000 --seed 11 --warmup 1000 --reps 10 --outdir runs/time/l100_m110_s11/
