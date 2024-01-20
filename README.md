# Thesis
## How to use
* The model has already been trained and saved as "modelpara.pt".Thus we don't need to train the model again.

  After modifying **input_string** in `main.py`, translate Germany into English using the command below.
  
  ```python=
  python3 main.py
  # The terminal will show =>
  # The input string is: "your input".
  # The output string is: "your input in English".
  ```
  
* If you want to try another structure of the model, you only need to modify the `main.py` to change the parameter of the model
  and then train the model again.
  
  Make sure to cancel the annotation to re-train the transformer model.
