# Thesis
## How to use
* The model has already been trained and save as "modelpara.pt".Thus we don't need to train th emodel again.
  After modify **input_string** in `main.py`,use the command below to translate the Germany into English.
  ```bash=
  python3 main.py
  # The terminal will show
  # The input string is : "your input".
  # The output string is : "your input in English".
  ```
* If you want to try another structure of the model, only need to modify the `main.py` to change the parameter of the model
  and then train the model again.
  Make sure to cancel the annotation so as to re-trained the transformer model.