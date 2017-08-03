## Enabling tab completion

Adapted instructions from [qiime2]: https://github.com/qiime2/q2cli/blob/master/README.md

### Bash

To enable tab completion in Bash, run the following command or add it to your `.bashrc`/`.bash_profile`:

```bash
source shogun-complete.sh
```

### ZSH

To enable tab completion in ZSH, run the following commands or add them to your `.zshrc`:

```bash
autoload bashcompinit && bashcompinit && source shogun-complete.sh
```
