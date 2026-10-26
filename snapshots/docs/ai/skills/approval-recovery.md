# Skill: approval-recovery

自動承認の容量不足・タイムアウトやsandbox起動失敗で解析が止まる際、未実行の操作を明示し、必要なら手動承認への切替を速やかに案内する。

## 障害の区別と停止条件

- 安全上の拒否、承認基盤の障害、解析・入力データのエラーを区別する。
- `Selected model is at capacity`や承認レビューtimeoutでは、操作は未実行の可能性がある。実行済み・保存済みと報告しない。
- `bwrap: Creating new namespace failed: No space left on device`はsandbox起動失敗。PicoDstやディスク容量の問題と断定しない。
- sandbox失敗時は正式なescalation経路で同じ限定操作の承認を要求する。レビュー障害は最大1回再試行できるが、復旧しなければ直ちにユーザーへ通知する。長い再試行ループを続けない。
- 複数agentで同じ障害への要求を繰り返さず、親agentへ集約して影響する処理を停止する。
- 明示的な安全上の拒否は別扱い。別コマンド・ツール等で同じ結果へ迂回しない。

## 手動承認への切替案内

未実行の操作、完了済みの操作、継続できる範囲を分けて知らせる。Todoを書けない場合もその旨を明示する。

VS CodeのCodex拡張では、入力欄下の権限メニューに表示があれば、ユーザー自身に`Approve for me`から`Ask for approval`へ変更してもらう。必要なら歯車→Codex Settings→Open config.tomlを案内する。

管理設定が許す場合の対応する値は以下。既存キーは重複追加せず置換し、ユーザー自身が変更する。

```toml
approval_policy = "on-request"
approvals_reviewer = "user"
sandbox_mode = "workspace-write"
```

agent自身が承認設定・管理ポリシーを変更しない。Full access・sandbox無効化・恒久的な広範囲許可を回避策として求めない。

## 手動中の実行と復帰

1. 読取りで途中出力・実際の変更・実行中処理を確認する。失敗した呼出しの部分実行も疑い、二重実行を防ぐ。
2. 当面必要な読取り・限定編集・移動・検証を整理して提示する。未確定の将来処理を包括承認へ含めない。
3. 正式な承認経路で実行し、結果を検証してTodo・ログへ記録する。
4. 当面の承認対象が終わったら、その旨を知らせ、自動承認へ戻すかをユーザーへ委ねる。
5. 既承認prefixの成功だけで、自動レビューやsandbox自体の復旧と断定しない。再発時はすぐ同じ案内を行う。

承認・サービス・ホスト資源の障害を完全に防止する保証はしない。目的は停滞を隠さず、早く安全な再開手順へ移ること。

## 根拠と保守

2026-10-19のユーザー指示に基づく。UI案内時は利用クライアントと公式資料を確認する。

- [Auto-review](https://developers.openai.com/codex/auto-review)
- [Sandbox・承認](https://learn.chatgpt.com/docs/sandboxing)
- [設定ファイル](https://developers.openai.com/codex/config-basic)

本手順の変更後は`script/sync_and_check_skills.sh`を実行する。
