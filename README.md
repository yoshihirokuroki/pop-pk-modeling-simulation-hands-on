# nlmixr2 & rxode2 母集団薬物動態解析ハンズオン

このリポジトリは、nlmixr2とrxode2パッケージを使用した母集団薬物動態（PopPK）解析のハンズオンセミナー用の教材です。

## 📚 概要

抗体医薬品（抗悪性腫瘍薬）を題材に、以下の内容を学びます：

-   データの可視化と探索的解析
-   1コンパートメントと2コンパートメントモデルの比較
-   共変量探索とモデル構築
-   モデル診断とバリデーション
-   （オプション）rxode2によるシミュレーション

## 🚀 クイックスタート（Posit Cloud）

### 方法1: GitHubから直接（推奨）

1.  [Posit Cloud](https://posit.cloud)にログイン
2.  **New Project** → **New Project from Git Repository**
3.  URL入力: <https://github.com/yoshihirokuroki/pop-pk-modeling-simulation-hands-on.git>
4.  プロジェクトが開いたら、Consoleで実行:

``` r
renv::restore()
```

> 注意⚠️ Posit CloudはLinuxで運用されており、renv::restore()でパッケージをインストールするためにソースコードからビルドするので非常に時間を要します。特にパッケージ{RcppEigen}はPosit Cloudに合理的な時間内でインストールできるかわかりません。実用的なハンズオン環境の目処が立つまでお待ちください。

### 方法2: テンプレートから

（テンプレートリンクを共有予定）

## 📁 ファイル構成

```         
nlmixr2-handson/
├── README.md                    # このファイル
├── .Rprofile                    # renv自動有効化
├── renv.lock                    # パッケージバージョン管理
├── presentation.qmd             # プレゼンテーション資料
├── handson.qmd                  # メインハンズオンファイル
├── data/
│   ├── pk_data.csv             # シミュレーションデータ
│   └── true_parameters.csv     # 真のパラメータ値（参考）
└── scripts/
    └── 01_generate_data.R      # データ生成スクリプト（参考）
```

## 🎯 学習目標

このハンズオンを完了すると、以下ができるようになります：

-   [ ] nlmixr2を使った基本的な母集団PK解析
-   [ ] モデル選択の実践（AIC、診断プロット）
-   [ ] 共変量効果の評価と解釈
-   [ ] 臨床的に意義のある結果の導出

## 📊 ケーススタディ

### 投与レジメン

-   **薬剤**: 抗体医薬品（架空）
-   **投与量**: 800 mg
-   **投与経路**: 1時間持続静脈内投与
-   **投与間隔**: 3週間ごと、6サイクル

### データセット

-   **被験者数**: 24名（男女比10:14）
-   **年齢範囲**: 35-85歳
-   **体重範囲**: 30-150 kg
-   **サンプリング**:
    -   Rich PK（サイクル1, 6）
    -   Sparse PK（サイクル2-5）

### 共変量

-   体重（WT）
-   年齢（AGE）
-   性別（SEX）
-   アルブミン（ALB）
-   ベースライン腫瘍径（TUMOR）

## 🔧 必要な環境

-   **R**: 4.3.2以上
-   **主要パッケージ**:
    -   nlmixr2 (2.1.2)
    -   rxode2 (2.1.2)
    -   dplyr, ggplot2, GGally, patchwork

すべてのパッケージは`renv::restore()`で自動インストールされます。

## 📖 使い方

### 1. 環境準備（5分）

Posit Cloudでプロジェクト作成後:

``` r
# パッケージインストール
renv::restore()

# データ確認
pk_data <- read.csv("data/pk_data.csv")
head(pk_data)
```

### 2. ハンズオン実施（40分）

`handson.qmd`を開き、各コードチャンクを実行しながら進めます。

**重要**: レンダリングはせず、**コードチャンクを個別に実行**してください。

### 3. オプション: rxode2シミュレーション（10分）

時間があれば、最終モデルを使ったシミュレーションを体験します。

## 🎓 前提知識

-   R言語の基本的な操作
-   薬物動態学の基礎知識
-   コンパートメントモデルの概念

## 📚 参考資料

### オンラインリソース

-   [nlmixr2公式サイト](https://nlmixr2.org/)
-   [rxode2ドキュメント](https://nlmixr2.github.io/rxode2/)
-   [Posit Cloud ガイド](https://posit.cloud/learn/guide)

### 推奨教科書

-   Gabrielsson & Weiner: *Pharmacokinetic and Pharmacodynamic Data Analysis*
-   Bonate: *Pharmacokinetic-Pharmacodynamic Modeling and Simulation*

## ❓ トラブルシューティング

### パッケージインストールエラー

``` r
# renvキャッシュをクリア
renv::purge()
renv::restore()
```

### フィッティングが収束しない

-   初期値を調整
-   `saemControl()`のパラメータを変更
-   アルゴリズムを変更（SAEM → FOCEI）

### メモリ不足

Posit Cloudの無料プランでは1GBのRAM制限があります。必要に応じて有料プランへのアップグレードを検討してください。

## 📞 サポート

質問や問題があれば：

-   Issueを作成
-   メール: [your.email\@example.com](mailto:your.email@example.com){.email}
-   セミナー中は講師に直接質問

## 📄 ライセンス

このプロジェクトはMITライセンスの下で公開されています。

## 🙏 謝辞

nlmixr2およびrxode2の開発チームに感謝します。

------------------------------------------------------------------------

**Happy Modeling!** 🎉

*最終更新: 2025年1月*
