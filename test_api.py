import time

import requests

BASE_URL = "http://localhost:8000/api"


def test_lifecycle():
    print("Testing Start...")
    r = requests.post(f"{BASE_URL}/control/start")
    print(f"Start response: {r.status_code}, {r.text}")

    time.sleep(2)

    print("Testing Status...")
    r = requests.get(f"{BASE_URL}/status")
    print(f"Status response: {r.status_code}, {r.text}")

    print("Testing Stop...")
    r = requests.post(f"{BASE_URL}/control/stop")
    print(f"Stop response: {r.status_code}, {r.text}")

    print("Testing Reset...")
    r = requests.post(f"{BASE_URL}/control/reset")
    print(f"Reset response: {r.status_code}, {r.text}")

    print("Testing Restart...")
    r = requests.post(f"{BASE_URL}/control/start")
    print(f"Restart response: {r.status_code}, {r.text}")

    time.sleep(1)
    print("Final Status check...")
    r = requests.get(f"{BASE_URL}/status")
    print(f"Final Status: {r.status_code}, {r.text}")


if __name__ == "__main__":
    test_lifecycle()
