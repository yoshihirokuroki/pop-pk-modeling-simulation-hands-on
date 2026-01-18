# ==============================================================================
# 抗体医薬品PK データ生成スクリプト
# ==============================================================================
# 目的: nlmixr2ハンズオン用のシミュレーションデータセット作成
# モデル: 2コンパートメント、比例誤差、M3法対応（BLQ）
# 共変量: アルブミン、性別、体重、年齢、ベースライン腫瘍径
# ==============================================================================

library(rxode2)
library(dplyr)
library(tidyr)

set.seed(123)

# ==============================================================================
# 1. 真のモデルパラメータ設定
# ==============================================================================

# 典型値（70kg男性、アルブミン4.0 g/dL）
# ★修正: 抗体医薬品らしい半減期になるようパラメータ調整
TVCL <- 0.20    # L/day (0.25 → 0.20: クリアランスを20%減少)
TVV1 <- 3.0     # L (中心コンパートメント、抗体は血管内に留まりやすい)
TVQ <- 0.15     # L/day (0.35 → 0.15: 分布を遅くして終末相を延長)
TVV2 <- 4.5     # L (2.8 → 4.5: 末梢容積を増やして貯蔵効果を高める)

# 個体間変動（対数正規分布、CV%）
# ★修正1: CLの個体間変動を小さくして共変量効果を検出しやすく
omega_CL <- 0.15  # 0.30 → 0.15 (共変量で説明される部分を除いた残差)
omega_V1 <- 0.25
omega_Q <- 0.35
omega_V2 <- 0.30

# 残差誤差（比例誤差）
sigma_prop <- 0.20

# 共変量効果
# CL = TVCL * (WT/70)^0.3 * (ALB/4.0)^(-1.2) * exp(eta_CL)
# ★修正2: アルブミン効果を大きくして検出しやすく
theta_WT_CL <- 0.3
theta_ALB_CL <- -1.2  # -0.6 → -1.2 (効果を2倍に)

# LLOQ設定
LLOQ <- 0.1  # μg/mL（抗体医薬品の高感度アッセイ）

# ==============================================================================
# 2. 被験者情報生成（24名）
# ==============================================================================

# n_subjects <- 16
n_subjects <- 24

# 性別（男女比55:45）
# sex <- rep(c("Male", "Female"), each = 8)
sex <- sample(
  x = c("Male", "Female"),
  size = n_subjects,
  replace = TRUE,
  prob = c(0.55, 0.45)
)

# 年齢（35-85歳）
age <- runif(n_subjects, 35, 85)

# 体重（性別と相関、男性：60-150kg、女性：30-90kg）
# 年齢との弱い負の相関も考慮
wt <- numeric(n_subjects)
for(i in 1:n_subjects) {
  if(sex[i] == "Male") {
    wt_mean <- 85 - (age[i] - 60) * 0.15  # 年齢による体重減少
    wt[i] <- rnorm(1, wt_mean, 15)
    wt[i] <- pmax(60, pmin(150, wt[i]))
  } else {
    wt_mean <- 60 - (age[i] - 60) * 0.12
    wt[i] <- rnorm(1, wt_mean, 12)
    wt[i] <- pmax(30, pmin(90, wt[i]))
  }
}

# アルブミン（2.8-4.5 g/dL）
# ベースライン腫瘍径との中程度の負の相関
# ★修正3: アルブミンの範囲を広げて共変量効果をより明確に
alb_base <- rnorm(n_subjects, 3.6, 0.6)  # 平均3.8→3.6、SD 0.4→0.6
alb_base <- pmax(2.5, pmin(4.8, alb_base))  # 範囲を拡大

# ベースライン腫瘍径（30-150 mm）
# アルブミンとの相関構造
tumor_size <- numeric(n_subjects)
for(i in 1:n_subjects) {
  # アルブミンが低い→腫瘍が大きい傾向
  tumor_mean <- 90 - (alb_base[i] - 3.8) * 25
  tumor_size[i] <- rnorm(1, tumor_mean, 20)
  tumor_size[i] <- pmax(30, pmin(150, tumor_size[i]))
}

# 被験者データフレーム
subjects <- data.frame(
  ID = 1:n_subjects,
  SEX = sex,
  AGE = round(age, 1),
  WT = round(wt, 1),
  ALB = round(alb_base, 2),
  TUMOR = round(tumor_size, 1)
)

# ==============================================================================
# 3. 投与レジメンとサンプリング時刻設定
# ==============================================================================

# 投与: 800 mg、1時間持続点滴、3週間（21日）ごと、6サイクル
dose_amt <- 800  # mg
infusion_duration <- 1  # hour
cycle_interval <- 21  # days

# サンプリングスケジュール
# Cycle 1 & 6: Rich sampling
rich_times <- c(0, 0.5, 1, 1.5, 2, 4, 8, 24, 72, 168, 336, 504)  # hours after dose

# Cycle 2-5: Sparse sampling
sparse_times <- c(0, 168, 504)  # Pre-dose, Day 7, Day 21

# 各被験者の投与・サンプリングデータ作成
create_subject_data <- function(id_val, subj_info) {
  
  dose_records <- data.frame(
    ID = id_val,
    TIME = (0:5) * cycle_interval * 24,  # 各サイクル開始時刻（hours）
    AMT = dose_amt,
    RATE = dose_amt / infusion_duration,
    DV = NA,
    EVID = 1,
    CMT = 1,
    # MDV = 1,
    BLQ = 0
  )
  
  # サンプリング
  sampling_records <- list()
  
  for(cycle in 1:6) {
    cycle_start <- (cycle - 1) * cycle_interval * 24
    
    if(cycle == 1 || cycle == 6) {
      sample_times <- rich_times
    } else {
      sample_times <- sparse_times
    }
    
    sampling_records[[cycle]] <- data.frame(
      ID = id_val,
      TIME = cycle_start + sample_times,
      AMT = 0,
      RATE = 0,
      DV = NA,
      EVID = 0,
      CMT = 2,  # 観測コンパートメント
      # MDV = 0,
      BLQ = 0
    )
  }
  
  sampling_df <- do.call(rbind, sampling_records)
  
  data <- rbind(dose_records, sampling_df)
  data <- data[order(data$TIME), ]
  
  # 共変量追加
  data$SEX <- subj_info$SEX
  data$AGE <- subj_info$AGE
  data$WT <- subj_info$WT
  data$ALB <- subj_info$ALB
  data$TUMOR <- subj_info$TUMOR
  
  return(data)
}

# 全被験者データ結合
all_data <- do.call(rbind, lapply(1:n_subjects, function(i) {
  create_subject_data(i, subjects[i, ])
}))

# ==============================================================================
# 4. rxode2モデル定義と濃度シミュレーション
# ==============================================================================

pk_model <- rxode2({
  # 共変量効果
  CL_cov = (WT/70)^theta_WT_CL * (ALB/4.0)^theta_ALB_CL
  
  # 個体パラメータ
  CL = TVCL * CL_cov * exp(eta_CL)
  V1 = TVV1 * exp(eta_V1)
  Q = TVQ * exp(eta_Q)
  V2 = TVV2 * exp(eta_V2)
  
  # 微分方程式
  d/dt(A1) = -CL/V1*A1 - Q/V1*A1 + Q/V2*A2
  d/dt(A2) = Q/V1*A1 - Q/V2*A2
  
  # 血中濃度（中心コンパートメント）
  Cp = A1/V1
  
  # 観測値（比例誤差）
  # Cp_obs = Cp * (1 + eps_prop)
})

# 個体間変動パラメータ生成
eta_matrix <- matrix(0, nrow = n_subjects, ncol = 4)
colnames(eta_matrix) <- c("eta_CL", "eta_V1", "eta_Q", "eta_V2")

for(i in 1:n_subjects) {
  eta_matrix[i, "eta_CL"] <- rnorm(1, 0, omega_CL)
  eta_matrix[i, "eta_V1"] <- rnorm(1, 0, omega_V1)
  eta_matrix[i, "eta_Q"] <- rnorm(1, 0, omega_Q)
  eta_matrix[i, "eta_V2"] <- rnorm(1, 0, omega_V2)
}

# シミュレーション実行
# params <- data.frame(
#   TVCL = TVCL,
#   TVV1 = TVV1,
#   TVQ = TVQ,
#   TVV2 = TVV2,
#   theta_WT_CL = theta_WT_CL,
#   theta_ALB_CL = theta_ALB_CL
# )
base_params <- c(
  TVCL = TVCL,
  TVV1 = TVV1,
  TVQ  = TVQ,
  TVV2 = TVV2,
  theta_WT_CL = theta_WT_CL,
  theta_ALB_CL = theta_ALB_CL
)

sim_data <- all_data
sim_results <- list()

for(i in 1:n_subjects) {
  subj_data <- sim_data[sim_data$ID == i, ]
  
  subj_params <- c(
    base_params,
    eta_matrix[i, ]   # ← これはすでに named numeric
  )

  events <- subj_data[, c("TIME", "AMT", "RATE", "EVID", "CMT", "WT", "ALB")]
  
  sim <- rxSolve(pk_model, 
                 params = subj_params,
                 events = events,
                 inits = c(A1 = 0, A2 = 0))
  
  sim_results[[i]] <- sim
}

# シミュレーション結果を元データに結合
for(i in 1:n_subjects) {
  subj_idx <- which(sim_data$ID == i & sim_data$EVID == 0)
  
  sim_conc <- sim_results[[i]]$Cp[sim_results[[i]]$time %in% 
                                    sim_data$TIME[subj_idx]]
  
  # 残差誤差追加
  eps <- rnorm(length(sim_conc), 0, sigma_prop)
  obs_conc <- sim_conc * (1 + eps)
  obs_conc <- pmax(0, obs_conc)  # 負の値を0に
  
  sim_data$DV[subj_idx] <- obs_conc
}

# ==============================================================================
# 5. BLQフラグ設定（M3法対応）
# ==============================================================================

sim_data$BLQ <- ifelse(sim_data$EVID == 0 & sim_data$DV < LLOQ, 1, 0)
sim_data$DV[sim_data$BLQ == 1] <- NA  # BLQはNAに
sim_data$LLOQ <- LLOQ
sim_data$LLOQ[sim_data$EVID == 1] <- NA_real_

# ==============================================================================
# 6. NONMEM形式への整形
# ==============================================================================

final_data <- sim_data %>%
  mutate(
    DV = round(DV, 3),
    SEXN = ifelse(SEX == "Male", 0, 1),  # 数値化
    CYCLE = ceiling(TIME / (cycle_interval * 24))
  ) %>%
  select(ID, TIME, AMT, RATE, DV, EVID, CMT, BLQ, LLOQ,
         WT, AGE, SEX, SEXN, ALB, TUMOR, CYCLE) %>%
  arrange(ID, TIME)

# ==============================================================================
# 7. データ保存
# ==============================================================================

# dataディレクトリ作成
if(!dir.exists("data")) dir.create("data")

write.csv(final_data, "data/pk_data.csv", row.names = FALSE, quote = FALSE)

# 真のパラメータ値も保存（参考用）
true_params <- data.frame(
  Parameter = c("TVCL", "TVV1", "TVQ", "TVV2", 
                "omega_CL", "omega_V1", "omega_Q", "omega_V2",
                "theta_WT_CL", "theta_ALB_CL", "sigma_prop", "LLOQ"),
  Value = c(TVCL, TVV1, TVQ, TVV2,
            omega_CL, omega_V1, omega_Q, omega_V2,
            theta_WT_CL, theta_ALB_CL, sigma_prop, LLOQ),
  Description = c(
    "Typical CL (L/day)", "Typical V1 (L)", "Typical Q (L/day)", "Typical V2 (L)",
    "IIV on CL (CV)", "IIV on V1 (CV)", "IIV on Q (CV)", "IIV on V2 (CV)",
    "WT exponent on CL", "ALB exponent on CL", "Proportional error (CV)", 
    "Lower limit of quantification (μg/mL)"
  )
)

write.csv(true_params, "data/true_parameters.csv", row.names = FALSE)

# 要約統計量表示
cat("\n=== データ生成完了 ===\n")
cat("被験者数:", n_subjects, "\n")
cat("総観測数:", sum(final_data$EVID == 0), "\n")
cat("BLQ数:", sum(final_data$BLQ == 1, na.rm = TRUE), "\n")
cat("BLQ割合:", 
    round(sum(final_data$BLQ == 1, na.rm = TRUE) / sum(final_data$EVID == 0) * 100, 1), "%\n")

print(summary(subjects[, c("AGE", "WT", "ALB", "TUMOR")]))
table(subjects$SEX)

