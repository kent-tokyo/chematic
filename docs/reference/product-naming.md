# Product naming

The official product names are `chematic` and `chematic-draw`. Both stay
lowercase, including at the beginning of a sentence. Do not add an extra
leading `s`, and do not capitalize the names as ordinary title-case words.

CI checks current public text and package metadata with:

```bash
python3 scripts/checks/check_product_names.py
```

Historical benchmark records, archived documents, internal identifiers, and
frozen wire-format values are outside this branding check. They remain stable
so naming cleanup cannot rewrite evidence or break compatibility contracts.
