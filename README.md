# Sochisirius_aug2025

В этом коде используются [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
Для создания и воспроизведения научного проекта под названием
> Sochisirius_aug2025

Для (локального) воспроизведения этого проекта:

1. Откройте Julia и выполните:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

При этом будут установлены все необходимые пакеты, чтобы вы могли запускать скрипты, и
все должно работать "из коробки", включая правильный поиск локальных путей.

Перед запуском скриптов выполните:
```julia
using DrWatson
@quickactivate "Sochisirius_aug2025"
```
которые автоматически активируют проект и включают обработку локального пути из Drwatson.

Непосредственно для запуска скрипта используйте:
```julia
include(scriptsdir("Название_скрипта.jl"))
```

При возникновении ошибки связанной с OPENSSL перезапустите скрипт не выходя из julia.

Для запуска многопоточной программы запускайте julia с опцией --threads=количество_потоков
