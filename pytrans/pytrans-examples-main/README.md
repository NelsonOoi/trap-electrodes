# Pytrans Models

Examples and showcase repo for waveform generation in actual traps with pytrans

## What is this repo?

This is not an active project, to generate waveforms for existing traps, even though some of the traps illustrated here are real. Its aim is

1) To provide working examples illustrating `pytrans` with both toy models and real use cases.
2) To build a library of documented models of traps that have used `pytrans`. Models should be uploaded here for documentation purpose, referencing their data source (or how to recreate it) and should refer to the respective project repositories for specific use cases.
2) To suggest a folder structure for future projects. Sort of a cookiecutter.

```
my-pytrans-project/
├─ deps/
│  ├─ pytrans/
│  ├─ .../
│
├─ model/
│  ├─ model.py
│  ├─ __init__.py
│
├─ scripts/
│  ├─ ...
├─ notebooks/
│  ├─ ...
│
├─ pyproject.toml
```

## Installation

git clone && `poetry install`